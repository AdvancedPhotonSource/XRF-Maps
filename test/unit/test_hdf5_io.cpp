// Unit tests for io::file::HDF5_IO and the free template helpers in hdf5_io.h
//
// Covers:
//   - parse_str_val_to_int<T>     (pure string-parsing logic, no I/O)
//   - translate_back_sens_num<T>  (lookup table, no I/O)
//   - translate_back_sens_unit<T> (lookup table, no I/O)
//   - HDF5_IO::load_spectra_volume           (error and success paths)
//   - HDF5_IO::load_spectra_volume_with_callback (error and success paths)
//
// Build with BUILD_TESTS=ON:
//   cmake -DBUILD_TESTS=ON ..
//   cmake --build . --target test_hdf5_io
//   ctest  (or ./bin/test_hdf5_io)

#include <gtest/gtest.h>
#include <hdf5.h>

#include <filesystem>
#include <vector>

#include "io/file/hdf5_io.h"
#include "data_struct/spectra_volume.h"
#include "data_struct/spectra.h"

using namespace io::file;

// ============================================================
// Helpers — create minimal MAPS_RAW-layout HDF5 files
// ============================================================

struct MapsRawParams
{
    size_t spectra_size  = 16;   // small to keep tests fast
    size_t rows          = 2;
    size_t cols          = 3;
    size_t n_detectors   = 4;
    size_t detector_num  = 0;
    float  livetime      = 0.10f;
    float  realtime      = 0.12f;
    float  input_counts  = 100.0f;
    float  output_counts = 100.0f;
    // When true, the spectra data[s,row,col] = float(s+1); otherwise uniform 5.0f
    bool   channel_ramp  = false;
};

static const char* detector_dset_name(size_t det)
{
    switch (det)
    {
        case 0: return "data_a";
        case 1: return "data_b";
        case 2: return "data_c";
        case 3: return "data_d";
        default: return "";
    }
}

// Returns 0 on success, negative on failure.
static int make_maps_raw_h5(const std::string& path, const MapsRawParams& p)
{
    hid_t fid = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (fid < 0) return -1;

    hid_t grp = H5Gcreate2(fid, "MAPS_RAW", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (grp < 0) { H5Fclose(fid); return -1; }

    // --- spectra dataset [spectra_size, rows, cols] ---
    {
        hsize_t dims[3] = { p.spectra_size, p.rows, p.cols };
        hid_t space = H5Screate_simple(3, dims, nullptr);
        const char* dname = detector_dset_name(p.detector_num);
        hid_t dset = H5Dcreate2(grp, dname, H5T_NATIVE_FLOAT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        size_t total = p.spectra_size * p.rows * p.cols;
        std::vector<float> buf(total);
        for (size_t s = 0; s < p.spectra_size; ++s)
            for (size_t r = 0; r < p.rows; ++r)
                for (size_t c = 0; c < p.cols; ++c)
                {
                    size_t idx = s * p.rows * p.cols + r * p.cols + c;
                    buf[idx] = p.channel_ramp ? float(s + 1) : 5.0f;
                }

        H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
        H5Dclose(dset);
        H5Sclose(space);
    }

    // --- metadata datasets [n_detectors, rows, cols] ---
    auto write_meta = [&](const char* name, float value) -> int
    {
        hsize_t dims[3] = { p.n_detectors, p.rows, p.cols };
        hid_t space = H5Screate_simple(3, dims, nullptr);
        hid_t dset  = H5Dcreate2(grp, name, H5T_NATIVE_FLOAT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        std::vector<float> buf(p.n_detectors * p.rows * p.cols, value);
        H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
        H5Dclose(dset);
        H5Sclose(space);
        return 0;
    };

    write_meta("livetime",     p.livetime);
    write_meta("realtime",     p.realtime);
    write_meta("inputcounts",  p.input_counts);
    write_meta("ouputcounts",  p.output_counts); // typo preserved from production code

    H5Gclose(grp);
    H5Fclose(fid);
    return 0;
}

// ============================================================
// Test fixture
// ============================================================

class HDF5_IO_Test : public ::testing::Test
{
protected:
    std::filesystem::path tmp_dir;

    void SetUp() override
    {
        tmp_dir = std::filesystem::temp_directory_path() / "xrf_hdf5_io_tests";
        std::filesystem::create_directories(tmp_dir);
    }

    void TearDown() override
    {
        std::filesystem::remove_all(tmp_dir);
    }

    std::string tmp(const std::string& name) const
    {
        return (tmp_dir / name).string();
    }

    static HDF5_IO* io() { return HDF5_IO::inst(); }
};

// ============================================================
// parse_str_val_to_int — free template function
// ============================================================

TEST(ParseStrValToInt, ReturnsNegOneWhenStartDelimiterMissing)
{
    EXPECT_EQ((parse_str_val_to_int<float>("START:", "\"", "no start here")), -1);
}

TEST(ParseStrValToInt, ReturnsNegOneWhenEndDelimiterMissing)
{
    // Start delimiter is found but end delimiter never appears
    EXPECT_EQ((parse_str_val_to_int<float>("\"", "\"", "\"42")), -1);
}

TEST(ParseStrValToInt, ParsesValueBetweenDelimiters)
{
    EXPECT_EQ((parse_str_val_to_int<float>("\"", "\"", "\"42\"rest")), 42);
}

TEST(ParseStrValToInt, ParsesValueWithCustomDelimiters)
{
    EXPECT_EQ((parse_str_val_to_int<float>("[", "]", "prefix[123]suffix")), 123);
}

TEST(ParseStrValToInt, ParsesZero)
{
    EXPECT_EQ((parse_str_val_to_int<float>("=", ";", "val=0;end")), 0);
}

// ============================================================
// translate_back_sens_num — free template function
// ============================================================

TEST(TranslateBackSensNum, AllValidInputsMappedCorrectly)
{
    EXPECT_FLOAT_EQ((translate_back_sens_num<float>(1)),   0.f);
    EXPECT_FLOAT_EQ((translate_back_sens_num<float>(2)),   1.f);
    EXPECT_FLOAT_EQ((translate_back_sens_num<float>(5)),   2.f);
    EXPECT_FLOAT_EQ((translate_back_sens_num<float>(10)),  3.f);
    EXPECT_FLOAT_EQ((translate_back_sens_num<float>(20)),  4.f);
    EXPECT_FLOAT_EQ((translate_back_sens_num<float>(50)),  5.f);
    EXPECT_FLOAT_EQ((translate_back_sens_num<float>(100)), 6.f);
    EXPECT_FLOAT_EQ((translate_back_sens_num<float>(200)), 7.f);
    EXPECT_FLOAT_EQ((translate_back_sens_num<float>(500)), 8.f);
}

TEST(TranslateBackSensNum, InvalidInputReturnsNegOne)
{
    EXPECT_FLOAT_EQ((translate_back_sens_num<float>(0)),   -1.f);
    EXPECT_FLOAT_EQ((translate_back_sens_num<float>(3)),   -1.f);
    EXPECT_FLOAT_EQ((translate_back_sens_num<float>(999)), -1.f);
}

TEST(TranslateBackSensNum, WorksWithDoubleType)
{
    EXPECT_DOUBLE_EQ((translate_back_sens_num<double>(100)), 6.0);
    EXPECT_DOUBLE_EQ((translate_back_sens_num<double>(7)),   -1.0);
}

// ============================================================
// translate_back_sens_unit — free template function
// ============================================================

TEST(TranslateBackSensUnit, AllValidUnitsMappedCorrectly)
{
    EXPECT_FLOAT_EQ((translate_back_sens_unit<float>("pA/V")), 0.f);
    EXPECT_FLOAT_EQ((translate_back_sens_unit<float>("nA/V")), 1.f);
    EXPECT_FLOAT_EQ((translate_back_sens_unit<float>("uA/V")), 2.f);
    EXPECT_FLOAT_EQ((translate_back_sens_unit<float>("mA/V")), 3.f);
}

TEST(TranslateBackSensUnit, UnknownUnitReturnsNegOne)
{
    EXPECT_FLOAT_EQ((translate_back_sens_unit<float>("")),        -1.f);
    EXPECT_FLOAT_EQ((translate_back_sens_unit<float>("kA/V")),    -1.f);
    EXPECT_FLOAT_EQ((translate_back_sens_unit<float>("unknown")), -1.f);
}

TEST(TranslateBackSensUnit, WorksWithDoubleType)
{
    EXPECT_DOUBLE_EQ((translate_back_sens_unit<double>("nA/V")),    1.0);
    EXPECT_DOUBLE_EQ((translate_back_sens_unit<double>("invalid")), -1.0);
}

// ============================================================
// load_spectra_volume — error paths
// ============================================================

TEST_F(HDF5_IO_Test, LoadSpectraVolume_NonExistentFile_ReturnsFalse)
{
    data_struct::Spectra_Volume<float> vol;
    vol.resize_and_zero(2, 3, 16);
    EXPECT_FALSE(io()->load_spectra_volume("/no/such/file.h5", 0, &vol));
}

TEST_F(HDF5_IO_Test, LoadSpectraVolume_MissingMapsRawGroup_ReturnsFalse)
{
    // Create a valid HDF5 file but without the /MAPS_RAW group
    std::string path = tmp("no_maps_raw.h5");
    hid_t fid = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    ASSERT_GE(fid, 0);
    H5Fclose(fid);

    data_struct::Spectra_Volume<float> vol;
    vol.resize_and_zero(2, 3, 16);
    EXPECT_FALSE(io()->load_spectra_volume(path, 0, &vol));
}

TEST_F(HDF5_IO_Test, LoadSpectraVolume_MissingDetectorDataset_ReturnsFalse)
{
    // /MAPS_RAW exists but the detector dataset (data_a) is absent
    std::string path = tmp("no_data_a.h5");
    hid_t fid = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    ASSERT_GE(fid, 0);
    hid_t grp = H5Gcreate2(fid, "MAPS_RAW", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(grp);
    H5Fclose(fid);

    data_struct::Spectra_Volume<float> vol;
    vol.resize_and_zero(2, 3, 16);
    EXPECT_FALSE(io()->load_spectra_volume(path, 0, &vol));
}

// ============================================================
// load_spectra_volume — success paths
// ============================================================

TEST_F(HDF5_IO_Test, LoadSpectraVolume_ValidFile_ReturnsTrue)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 0;
    std::string path = tmp("maps_raw_ok.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    data_struct::Spectra_Volume<float> vol;
    vol.resize_and_zero(p.rows, p.cols, p.spectra_size);
    EXPECT_TRUE(io()->load_spectra_volume(path, 0, &vol));
}

TEST_F(HDF5_IO_Test, LoadSpectraVolume_ValidFile_CorrectDimensions)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 0;
    std::string path = tmp("maps_raw_dims.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    data_struct::Spectra_Volume<float> vol;
    vol.resize_and_zero(p.rows, p.cols, p.spectra_size);
    ASSERT_TRUE(io()->load_spectra_volume(path, 0, &vol));

    EXPECT_EQ(vol.rows(), p.rows);
    EXPECT_EQ(vol.cols(), p.cols);
    EXPECT_EQ(vol.samples_size(), p.spectra_size);
}

TEST_F(HDF5_IO_Test, LoadSpectraVolume_UniformFill_CorrectSpectraValues)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 0;
    p.channel_ramp = false;  // all channels = 5.0f
    std::string path = tmp("maps_raw_uniform.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    data_struct::Spectra_Volume<float> vol;
    vol.resize_and_zero(p.rows, p.cols, p.spectra_size);
    ASSERT_TRUE(io()->load_spectra_volume(path, 0, &vol));

    for (size_t row = 0; row < p.rows; ++row)
        for (size_t col = 0; col < p.cols; ++col)
            for (size_t s = 0; s < p.spectra_size; ++s)
                EXPECT_FLOAT_EQ(vol[row][col][s], 5.0f)
                    << "row=" << row << " col=" << col << " s=" << s;
}

TEST_F(HDF5_IO_Test, LoadSpectraVolume_ChannelRamp_CorrectSpectraValues)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 0;
    p.channel_ramp = true;  // channel s has value float(s+1)
    std::string path = tmp("maps_raw_ramp.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    data_struct::Spectra_Volume<float> vol;
    vol.resize_and_zero(p.rows, p.cols, p.spectra_size);
    ASSERT_TRUE(io()->load_spectra_volume(path, 0, &vol));

    for (size_t row = 0; row < p.rows; ++row)
        for (size_t col = 0; col < p.cols; ++col)
            for (size_t s = 0; s < p.spectra_size; ++s)
                EXPECT_FLOAT_EQ(vol[row][col][s], float(s + 1))
                    << "row=" << row << " col=" << col << " s=" << s;
}

TEST_F(HDF5_IO_Test, LoadSpectraVolume_ValidFile_CorrectRealtime)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 0;
    p.realtime = 0.12f;
    p.input_counts = 100.0f;  p.output_counts = 100.0f;  // equal → livetime == realtime
    std::string path = tmp("maps_raw_timing.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    data_struct::Spectra_Volume<float> vol;
    vol.resize_and_zero(p.rows, p.cols, p.spectra_size);
    ASSERT_TRUE(io()->load_spectra_volume(path, 0, &vol));

    for (size_t row = 0; row < p.rows; ++row)
        for (size_t col = 0; col < p.cols; ++col)
            EXPECT_NEAR(vol[row][col].elapsed_realtime(), p.realtime, 1e-5f)
                << "row=" << row << " col=" << col;
}

TEST_F(HDF5_IO_Test, LoadSpectraVolume_EqualInOut_LivetimeEqualsRealtime)
{
    // recalc_elapsed_livetime: lt = rt * out/in; equal in/out → lt == rt
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 0;
    p.realtime = 0.12f;  p.livetime = 0.10f;
    p.input_counts = 200.0f;  p.output_counts = 200.0f;
    std::string path = tmp("maps_raw_lt_recalc.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    data_struct::Spectra_Volume<float> vol;
    vol.resize_and_zero(p.rows, p.cols, p.spectra_size);
    ASSERT_TRUE(io()->load_spectra_volume(path, 0, &vol));

    for (size_t row = 0; row < p.rows; ++row)
        for (size_t col = 0; col < p.cols; ++col)
            // After recalc: livetime = realtime * out/in = 0.12 * 1.0 = 0.12
            EXPECT_NEAR(vol[row][col].elapsed_livetime(), p.realtime, 1e-5f)
                << "row=" << row << " col=" << col;
}

TEST_F(HDF5_IO_Test, LoadSpectraVolume_Detector1_LoadsData)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 1;
    p.channel_ramp = true;
    std::string path = tmp("maps_raw_det1.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    data_struct::Spectra_Volume<float> vol;
    vol.resize_and_zero(p.rows, p.cols, p.spectra_size);
    ASSERT_TRUE(io()->load_spectra_volume(path, 1, &vol));

    for (size_t s = 0; s < p.spectra_size; ++s)
        EXPECT_FLOAT_EQ(vol[0][0][s], float(s + 1)) << "s=" << s;
}

// ============================================================
// load_spectra_volume_with_callback — error paths
// ============================================================

TEST_F(HDF5_IO_Test, LoadSpectraVolumeWithCallback_NonExistentFile_ReturnsFalse)
{
    data_struct::IO_Callback_Func_Def<float> cb =
        [](size_t, size_t, size_t, size_t, size_t, data_struct::Spectra<float>*, void*) {};

    EXPECT_FALSE(io()->load_spectra_volume_with_callback("/no/such/file.h5", {0}, cb, nullptr));
}

TEST_F(HDF5_IO_Test, LoadSpectraVolumeWithCallback_MissingMapsRawGroup_ReturnsFalse)
{
    std::string path = tmp("cb_no_maps_raw.h5");
    hid_t fid = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    ASSERT_GE(fid, 0);
    H5Fclose(fid);

    data_struct::IO_Callback_Func_Def<float> cb =
        [](size_t, size_t, size_t, size_t, size_t, data_struct::Spectra<float>*, void*) {};

    EXPECT_FALSE(io()->load_spectra_volume_with_callback(path, {0}, cb, nullptr));
}

// ============================================================
// load_spectra_volume_with_callback — success paths
// ============================================================

TEST_F(HDF5_IO_Test, LoadSpectraVolumeWithCallback_ValidFile_ReturnsTrue)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 0;
    std::string path = tmp("cb_ok.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    data_struct::IO_Callback_Func_Def<float> cb =
        [](size_t, size_t, size_t, size_t, size_t, data_struct::Spectra<float>* s, void*)
        { delete s; };

    EXPECT_TRUE(io()->load_spectra_volume_with_callback(path, {0}, cb, nullptr));
}

TEST_F(HDF5_IO_Test, LoadSpectraVolumeWithCallback_CallbackInvokedRowsTimesCols)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 3;  p.cols = 4;  p.detector_num = 0;
    std::string path = tmp("cb_count.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    size_t call_count = 0;
    data_struct::IO_Callback_Func_Def<float> cb =
        [&](size_t, size_t, size_t, size_t, size_t, data_struct::Spectra<float>* s, void*)
        { ++call_count; delete s; };

    ASSERT_TRUE(io()->load_spectra_volume_with_callback(path, {0}, cb, nullptr));
    EXPECT_EQ(call_count, p.rows * p.cols);
}

TEST_F(HDF5_IO_Test, LoadSpectraVolumeWithCallback_MultipleDetectors_CallbackCountScales)
{
    // Create two detector datasets in one file and request both.
    const size_t rows = 2, cols = 3, spectra_size = 16;

    std::string path = tmp("cb_multi_det.h5");
    hid_t fid = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    ASSERT_GE(fid, 0);
    hid_t grp = H5Gcreate2(fid, "MAPS_RAW", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Write data_a and data_b
    {
        hsize_t dims[3] = { spectra_size, rows, cols };
        hid_t space = H5Screate_simple(3, dims, nullptr);
        std::vector<float> buf(spectra_size * rows * cols, 1.0f);
        for (const char* nm : {"data_a", "data_b"})
        {
            hid_t d = H5Dcreate2(grp, nm, H5T_NATIVE_FLOAT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(d, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
            H5Dclose(d);
        }
        H5Sclose(space);
    }

    // Write metadata [4, rows, cols]
    {
        hsize_t dims[3] = { 4, rows, cols };
        hid_t space = H5Screate_simple(3, dims, nullptr);
        std::vector<float> buf(4 * rows * cols, 0.1f);
        for (const char* nm : {"livetime", "realtime", "inputcounts", "ouputcounts"})
        {
            hid_t d = H5Dcreate2(grp, nm, H5T_NATIVE_FLOAT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(d, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
            H5Dclose(d);
        }
        H5Sclose(space);
    }

    H5Gclose(grp);
    H5Fclose(fid);

    size_t call_count = 0;
    data_struct::IO_Callback_Func_Def<float> cb =
        [&](size_t, size_t, size_t, size_t, size_t, data_struct::Spectra<float>* s, void*)
        { ++call_count; delete s; };

    ASSERT_TRUE(io()->load_spectra_volume_with_callback(path, {0, 1}, cb, nullptr));
    // rows * cols * n_detectors = 2 * 3 * 2 = 12
    EXPECT_EQ(call_count, rows * cols * 2u);
}

TEST_F(HDF5_IO_Test, LoadSpectraVolumeWithCallback_CorrectDetectorIdInCallback)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 0;
    std::string path = tmp("cb_det_id.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    std::vector<size_t> det_ids;
    data_struct::IO_Callback_Func_Def<float> cb =
        [&](size_t, size_t, size_t, size_t, size_t det, data_struct::Spectra<float>* s, void*)
        { det_ids.push_back(det); delete s; };

    ASSERT_TRUE(io()->load_spectra_volume_with_callback(path, {0}, cb, nullptr));
    ASSERT_FALSE(det_ids.empty());
    for (size_t id : det_ids)
        EXPECT_EQ(id, 0u);
}

TEST_F(HDF5_IO_Test, LoadSpectraVolumeWithCallback_CorrectRowColIndices)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 0;
    std::string path = tmp("cb_rowcol.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    std::vector<std::pair<size_t, size_t>> row_cols;
    data_struct::IO_Callback_Func_Def<float> cb =
        [&](size_t row, size_t col, size_t, size_t, size_t, data_struct::Spectra<float>* s, void*)
        { row_cols.emplace_back(row, col); delete s; };

    ASSERT_TRUE(io()->load_spectra_volume_with_callback(path, {0}, cb, nullptr));
    ASSERT_EQ(row_cols.size(), p.rows * p.cols);

    size_t idx = 0;
    for (size_t r = 0; r < p.rows; ++r)
        for (size_t c = 0; c < p.cols; ++c, ++idx)
        {
            EXPECT_EQ(row_cols[idx].first,  r) << "idx=" << idx;
            EXPECT_EQ(row_cols[idx].second, c) << "idx=" << idx;
        }
}

TEST_F(HDF5_IO_Test, LoadSpectraVolumeWithCallback_CorrectSpectraDataInCallback)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 0;
    p.channel_ramp = true;
    std::string path = tmp("cb_spectra.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    std::vector<data_struct::Spectra<float>> received;
    data_struct::IO_Callback_Func_Def<float> cb =
        [&](size_t, size_t, size_t, size_t, size_t, data_struct::Spectra<float>* s, void*)
        {
            if (s) { received.push_back(*s); delete s; }
        };

    ASSERT_TRUE(io()->load_spectra_volume_with_callback(path, {0}, cb, nullptr));
    ASSERT_FALSE(received.empty());

    for (auto& spec : received)
        for (size_t s = 0; s < p.spectra_size; ++s)
            EXPECT_FLOAT_EQ(spec[s], float(s + 1)) << "s=" << s;
}

TEST_F(HDF5_IO_Test, LoadSpectraVolumeWithCallback_UserDataPassedThrough)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 1;  p.cols = 1;  p.detector_num = 0;
    std::string path = tmp("cb_userdata.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    int sentinel = 42;
    int* received_ptr = nullptr;
    data_struct::IO_Callback_Func_Def<float> cb =
        [&](size_t, size_t, size_t, size_t, size_t, data_struct::Spectra<float>* s, void* ud)
        { received_ptr = static_cast<int*>(ud); delete s; };

    ASSERT_TRUE(io()->load_spectra_volume_with_callback(path, {0}, cb, &sentinel));
    ASSERT_NE(received_ptr, nullptr);
    EXPECT_EQ(*received_ptr, 42);
}

TEST_F(HDF5_IO_Test, LoadSpectraVolumeWithCallback_TotalDimensionsPassedCorrectly)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 3;  p.cols = 5;  p.detector_num = 0;
    std::string path = tmp("cb_total_dims.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    size_t obs_total_rows = 0, obs_total_cols = 0;
    data_struct::IO_Callback_Func_Def<float> cb =
        [&](size_t, size_t, size_t tr, size_t tc, size_t, data_struct::Spectra<float>* s, void*)
        { obs_total_rows = tr;  obs_total_cols = tc;  delete s; };

    ASSERT_TRUE(io()->load_spectra_volume_with_callback(path, {0}, cb, nullptr));
    EXPECT_EQ(obs_total_rows, p.rows);
    EXPECT_EQ(obs_total_cols, p.cols);
}

TEST_F(HDF5_IO_Test, LoadSpectraVolumeWithCallback_TimingMetadataCorrect)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 1;  p.cols = 1;  p.detector_num = 0;
    p.realtime = 0.25f;  p.input_counts = 500.0f;  p.output_counts = 500.0f;
    std::string path = tmp("cb_timing.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    float obs_rt = -1.f;
    data_struct::IO_Callback_Func_Def<float> cb =
        [&](size_t, size_t, size_t, size_t, size_t, data_struct::Spectra<float>* s, void*)
        { if (s) { obs_rt = s->elapsed_realtime(); delete s; } };

    ASSERT_TRUE(io()->load_spectra_volume_with_callback(path, {0}, cb, nullptr));
    EXPECT_NEAR(obs_rt, p.realtime, 1e-5f);
}

// ============================================================
// load_spectra_volume — remaining detectors (2 and 3)
// ============================================================

TEST_F(HDF5_IO_Test, LoadSpectraVolume_Detector2_LoadsData)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 2;
    p.channel_ramp = true;
    std::string path = tmp("maps_raw_det2.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    data_struct::Spectra_Volume<float> vol;
    vol.resize_and_zero(p.rows, p.cols, p.spectra_size);
    ASSERT_TRUE(io()->load_spectra_volume(path, 2, &vol));

    for (size_t s = 0; s < p.spectra_size; ++s)
        EXPECT_FLOAT_EQ(vol[0][0][s], float(s + 1)) << "s=" << s;
}

TEST_F(HDF5_IO_Test, LoadSpectraVolume_Detector3_LoadsData)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 3;
    p.channel_ramp = true;
    std::string path = tmp("maps_raw_det3.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    data_struct::Spectra_Volume<float> vol;
    vol.resize_and_zero(p.rows, p.cols, p.spectra_size);
    ASSERT_TRUE(io()->load_spectra_volume(path, 3, &vol));

    for (size_t s = 0; s < p.spectra_size; ++s)
        EXPECT_FLOAT_EQ(vol[0][0][s], float(s + 1)) << "s=" << s;
}

// ============================================================
// load_spectra_volume — double precision type
// ============================================================

TEST_F(HDF5_IO_Test, LoadSpectraVolume_Double_ChannelRamp)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 0;
    p.channel_ramp = true;
    std::string path = tmp("maps_raw_double.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    data_struct::Spectra_Volume<double> vol;
    vol.resize_and_zero(p.rows, p.cols, p.spectra_size);
    ASSERT_TRUE(io()->load_spectra_volume(path, 0, &vol));

    for (size_t s = 0; s < p.spectra_size; ++s)
        EXPECT_NEAR(vol[0][0][s], double(s + 1), 1e-6) << "s=" << s;
}

TEST_F(HDF5_IO_Test, LoadSpectraVolume_Double_NonExistentFile_ReturnsFalse)
{
    data_struct::Spectra_Volume<double> vol;
    vol.resize_and_zero(2, 3, 16);
    EXPECT_FALSE(io()->load_spectra_volume<double>("/no/such/file.h5", 0, &vol));
}

// ============================================================
// load_spectra_volume_with_callback — double precision type
// ============================================================

TEST_F(HDF5_IO_Test, LoadSpectraVolumeWithCallback_Double_ChannelRamp)
{
    MapsRawParams p;
    p.spectra_size = 16;  p.rows = 2;  p.cols = 3;  p.detector_num = 0;
    p.channel_ramp = true;
    std::string path = tmp("cb_double.h5");
    ASSERT_EQ(0, make_maps_raw_h5(path, p));

    std::vector<data_struct::Spectra<double>> received;
    data_struct::IO_Callback_Func_Def<double> cb =
        [&](size_t, size_t, size_t, size_t, size_t, data_struct::Spectra<double>* s, void*)
        { if (s) { received.push_back(*s); delete s; } };

    ASSERT_TRUE(io()->load_spectra_volume_with_callback(path, {0}, cb, nullptr));
    ASSERT_FALSE(received.empty());

    for (auto& spec : received)
        for (size_t s = 0; s < p.spectra_size; ++s)
            EXPECT_NEAR(spec[s], double(s + 1), 1e-6) << "s=" << s;
}

// ============================================================
// parse_str_val_to_int — additional edge cases
// ============================================================

TEST(ParseStrValToInt, ParsesNegativeValue)
{
    EXPECT_EQ((parse_str_val_to_int<float>("[", "]", "val[-7]end")), -7);
}

TEST(ParseStrValToInt, ParsesMultiDigitValue)
{
    EXPECT_EQ((parse_str_val_to_int<float>("<", ">", "size<2048>px")), 2048);
}

TEST(ParseStrValToInt, StartDelimiterAtEndOfString_ReturnsFalse)
{
    // Start delimiter found but nothing after it
    EXPECT_EQ((parse_str_val_to_int<float>("[", "]", "value[")), -1);
}

// ============================================================
// load_spectra_volume_emd_with_callback — error paths
// ============================================================

TEST_F(HDF5_IO_Test, LoadSpectraVolumeEmdWithCallback_NonExistentFile_ReturnsFalse)
{
    data_struct::IO_Callback_Func_Def<float> cb =
        [](size_t, size_t, size_t, size_t, size_t, data_struct::Spectra<float>* s, void*) { delete s; };

    EXPECT_FALSE(io()->load_spectra_volume_emd_with_callback("/no/such/file.emd", {0}, cb, nullptr));
}

TEST_F(HDF5_IO_Test, LoadSpectraVolumeEmdWithCallback_EmptyH5File_ReturnsFalse)
{
    std::string path = tmp("empty_emd.h5");
    hid_t fid = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    ASSERT_GE(fid, 0);
    H5Fclose(fid);

    data_struct::IO_Callback_Func_Def<float> cb =
        [](size_t, size_t, size_t, size_t, size_t, data_struct::Spectra<float>* s, void*) { delete s; };

    EXPECT_FALSE(io()->load_spectra_volume_emd_with_callback(path, {0}, cb, nullptr));
}
