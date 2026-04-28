// Unit tests for functions in io/file/hl_file_io.h
//
// Covers:
//   - cb_load_spectra_data_helper  (pure in-memory accumulation logic)
//   - generate_fit_routine         (factory returning typed fit routines)
//   - load_override_params         (file-based parameter loading, error paths)
//   - load_and_integrate_spectra_volume  (MCA file path)
//
// Build with BUILD_TESTS=ON:
//   cmake -DBUILD_TESTS=ON ..
//   cmake --build . --target test_hl_file_io
//   ctest  (or ./bin/test_hl_file_io)

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "io/file/hl_file_io.h"
#include "data_struct/spectra.h"
#include "data_struct/params_override.h"
#include "fitting/routines/roi_fit_routine.h"
#include "fitting/routines/svd_fit_routine.h"
#include "fitting/routines/nnls_fit_routine.h"
#include "fitting/routines/param_optimized_fit_routine.h"
#include "fitting/routines/matrix_optimized_fit_routine.h"

namespace fs = std::filesystem;

using namespace data_struct;
using namespace fitting::routines;

// ============================================================
// cb_load_spectra_data_helper — pure accumulation logic, no I/O
// ============================================================

TEST(CbLoadSpectraDataHelper, NullUserData_DeletesSpectraWithoutCrash)
{
    // spectra is heap-allocated; helper must delete it even when user_data is null
    auto* spectra = new Spectra<double>(8);
    spectra->setConstant(1.0);
    io::file::cb_load_spectra_data_helper<double>(0, 0, 1, 1, 0, spectra, nullptr);
    // reaching here (and ASAN not firing) confirms the delete happened cleanly
}

TEST(CbLoadSpectraDataHelper, NullSpectra_DoesNotModifyIntegrated)
{
    Spectra<double> integrated(8);
    integrated.setConstant(3.0);
    io::file::cb_load_spectra_data_helper<double>(0, 0, 1, 1, 0, nullptr, &integrated);
    EXPECT_EQ(integrated.size(), 8);
    for (int i = 0; i < integrated.size(); ++i)
        EXPECT_DOUBLE_EQ(integrated[i], 3.0);
}

TEST(CbLoadSpectraDataHelper, BothNull_DoesNotCrash)
{
    io::file::cb_load_spectra_data_helper<double>(0, 0, 1, 1, 0, nullptr, nullptr);
}

TEST(CbLoadSpectraDataHelper, SameSize_AccumulatesIntoIntegrated)
{
    Spectra<double> integrated(4);
    integrated.setConstant(10.0);
    auto* spectra = new Spectra<double>(4);
    spectra->setConstant(5.0);

    io::file::cb_load_spectra_data_helper<double>(0, 0, 1, 1, 0, spectra, &integrated);

    ASSERT_EQ(integrated.size(), 4);
    for (int i = 0; i < integrated.size(); ++i)
        EXPECT_DOUBLE_EQ(integrated[i], 15.0);
}

TEST(CbLoadSpectraDataHelper, DifferentSize_CopiesSpectraIntoIntegrated)
{
    Spectra<double> integrated(4);
    integrated.setConstant(999.0);
    auto* spectra = new Spectra<double>(8);
    spectra->setConstant(3.0);

    io::file::cb_load_spectra_data_helper<double>(0, 0, 1, 1, 0, spectra, &integrated);

    ASSERT_EQ(integrated.size(), 8);
    for (int i = 0; i < integrated.size(); ++i)
        EXPECT_DOUBLE_EQ(integrated[i], 3.0);
}

TEST(CbLoadSpectraDataHelper, EmptyIntegrated_CopiesFirstSpectra)
{
    Spectra<double> integrated; // size 0
    auto* spectra = new Spectra<double>(4);
    spectra->setConstant(7.0);

    io::file::cb_load_spectra_data_helper<double>(0, 0, 1, 1, 0, spectra, &integrated);

    ASSERT_EQ(integrated.size(), 4);
    for (int i = 0; i < integrated.size(); ++i)
        EXPECT_DOUBLE_EQ(integrated[i], 7.0);
}

TEST(CbLoadSpectraDataHelper, MultipleAccumulationCalls_SumsCorrectly)
{
    Spectra<double> integrated(4);
    integrated.setZero();

    for (int call = 1; call <= 3; ++call)
    {
        auto* spectra = new Spectra<double>(4);
        spectra->setConstant(static_cast<double>(call));
        io::file::cb_load_spectra_data_helper<double>(0, 0, 1, 1, 0, spectra, &integrated);
    }

    // 1 + 2 + 3 = 6
    ASSERT_EQ(integrated.size(), 4);
    for (int i = 0; i < integrated.size(); ++i)
        EXPECT_DOUBLE_EQ(integrated[i], 6.0);
}

// ============================================================
// generate_fit_routine — factory, no I/O or optimizer state needed
// ============================================================

TEST(GenerateFitRoutine, ROI_ReturnsROIFitRoutine)
{
    auto* r = io::file::generate_fit_routine<double>(Fitting_Routines::ROI, nullptr);
    ASSERT_NE(r, nullptr);
    EXPECT_NE(dynamic_cast<ROI_Fit_Routine<double>*>(r), nullptr);
    delete r;
}

TEST(GenerateFitRoutine, SVD_ReturnsSVDFitRoutine)
{
    auto* r = io::file::generate_fit_routine<double>(Fitting_Routines::SVD, nullptr);
    ASSERT_NE(r, nullptr);
    EXPECT_NE(dynamic_cast<SVD_Fit_Routine<double>*>(r), nullptr);
    delete r;
}

TEST(GenerateFitRoutine, NNLS_ReturnsNNLSFitRoutine)
{
    auto* r = io::file::generate_fit_routine<double>(Fitting_Routines::NNLS, nullptr);
    ASSERT_NE(r, nullptr);
    EXPECT_NE(dynamic_cast<NNLS_Fit_Routine<double>*>(r), nullptr);
    delete r;
}

TEST(GenerateFitRoutine, GaussTails_ReturnsParamOptimizedNotMatrix)
{
    auto* r = io::file::generate_fit_routine<double>(Fitting_Routines::GAUSS_TAILS, nullptr);
    ASSERT_NE(r, nullptr);
    // Must be Param_Optimized, but NOT the Matrix subclass
    EXPECT_NE(dynamic_cast<Param_Optimized_Fit_Routine<double>*>(r), nullptr);
    EXPECT_EQ(dynamic_cast<Matrix_Optimized_Fit_Routine<double>*>(r), nullptr);
    delete r;
}

TEST(GenerateFitRoutine, GaussMatrix_ReturnsMatrixOptimizedFitRoutine)
{
    auto* r = io::file::generate_fit_routine<double>(Fitting_Routines::GAUSS_MATRIX, nullptr);
    ASSERT_NE(r, nullptr);
    EXPECT_NE(dynamic_cast<Matrix_Optimized_Fit_Routine<double>*>(r), nullptr);
    delete r;
}

TEST(GenerateFitRoutine, FloatInstantiation_ROI_ReturnsNonNull)
{
    auto* r = io::file::generate_fit_routine<float>(Fitting_Routines::ROI, nullptr);
    ASSERT_NE(r, nullptr);
    EXPECT_NE(dynamic_cast<ROI_Fit_Routine<float>*>(r), nullptr);
    delete r;
}

// ============================================================
// load_override_params — error paths only (no element DB loaded)
// ============================================================

TEST(LoadOverrideParams, NonexistentDirectory_ReturnsFalse)
{
    Params_Override<double> po;
    bool ok = io::file::load_override_params<double>("/nonexistent_xrf_test_dir_xyz/", -1, po);
    EXPECT_FALSE(ok);
}

TEST(LoadOverrideParams, EmptyDirectoryString_ReturnsFalse)
{
    Params_Override<double> po;
    bool ok = io::file::load_override_params<double>("", -1, po);
    EXPECT_FALSE(ok);
}

TEST(LoadOverrideParams, NonexistentDetectorFile_ReturnsFalse)
{
    Params_Override<double> po;
    bool ok = io::file::load_override_params<double>("/nonexistent_xrf_test_dir_xyz/", 0, po);
    EXPECT_FALSE(ok);
}

// ============================================================
// load_and_integrate_spectra_volume — MCA file path
// ============================================================

static const std::string s_tmp_dir = "/tmp/xrf_hl_file_io_tests";

// Writes a minimal version-3.1 MCA text file and returns its path.
static std::string make_minimal_mca(const std::string& dir,
                                    const std::string& name,
                                    int channels = 8)
{
    std::string path = dir + "/" + name;
    std::ofstream f(path);
    f << "VERSION: 3.1\n";
    f << "ELEMENTS: 1\n";
    f << "DATE: \n";
    f << "CHANNELS: " << channels << "\n";
    f << "REAL_TIME: 1.5\n";
    f << "LIVE_TIME: 1.2\n";
    f << "CAL_OFFSET: 0.0\n";
    f << "CAL_SLOPE: 1.0\n";
    f << "CAL_QUAD: 0.0\n";
    f << "DATA: \n";
    for (int i = 0; i < channels; ++i)
        f << static_cast<float>(i + 1) << "\n";
    return path;
}

class HlFileIoMcaTest : public ::testing::Test
{
protected:
    static void SetUpTestSuite()
    {
        fs::create_directories(s_tmp_dir);
    }

    static void TearDownTestSuite()
    {
        fs::remove_all(s_tmp_dir);
    }
};

TEST_F(HlFileIoMcaTest, LoadAndIntegrate_ValidMcaFile_ReturnsTrue)
{
    make_minimal_mca(s_tmp_dir, "scan_test.mca", 16);

    Spectra<double> integrated;
    Params_Override<double> po;

    bool ok = io::file::load_and_integrate_spectra_volume<double>(
        s_tmp_dir, "scan_test.mca", 0, integrated, &po);

    EXPECT_TRUE(ok);
    EXPECT_EQ(integrated.size(), 16);
    for (int i = 0; i < 16; ++i)
        EXPECT_FLOAT_EQ(static_cast<float>(integrated[i]), static_cast<float>(i + 1));
}

TEST_F(HlFileIoMcaTest, LoadAndIntegrate_MissingMcaFile_ReturnsFalse)
{
    Spectra<double> integrated;
    Params_Override<double> po;

    bool ok = io::file::load_and_integrate_spectra_volume<double>(
        s_tmp_dir, "does_not_exist.mca", 0, integrated, &po);

    EXPECT_FALSE(ok);
}

TEST_F(HlFileIoMcaTest, LoadAndIntegrate_Mca1Suffix_ReturnsTrue)
{
    make_minimal_mca(s_tmp_dir, "scan2.mca1", 8);

    Spectra<double> integrated;
    Params_Override<double> po;

    bool ok = io::file::load_and_integrate_spectra_volume<double>(
        s_tmp_dir, "scan2.mca1", 1, integrated, &po);

    EXPECT_TRUE(ok);
    EXPECT_EQ(integrated.size(), 8);
}

TEST_F(HlFileIoMcaTest, LoadAndIntegrate_SetsElapsedTimes)
{
    make_minimal_mca(s_tmp_dir, "scan_times.mca", 4);

    Spectra<double> integrated;
    Params_Override<double> po;

    io::file::load_and_integrate_spectra_volume<double>(
        s_tmp_dir, "scan_times.mca", 0, integrated, &po);

    EXPECT_FLOAT_EQ(static_cast<float>(integrated.elapsed_realtime()), 1.5f);
    EXPECT_FLOAT_EQ(static_cast<float>(integrated.elapsed_livetime()), 1.2f);
}
