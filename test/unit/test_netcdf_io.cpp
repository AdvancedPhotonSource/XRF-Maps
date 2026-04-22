// Unit tests for io::file::NetCDF_IO
//
// The reader expects a very specific binary layout packed into a NetCDF variable
// named "array_data".  The helpers below construct minimal files that match that
// layout so each public method can be exercised in isolation.
//
// Build with BUILD_TESTS=ON:
//   cmake -DBUILD_TESTS=ON ..
//   cmake --build . --target test_netcdf_io
//   ctest  (or ./bin/test_netcdf_io)

#include <gtest/gtest.h>
#include <netcdf.h>

#include <filesystem>
#include <functional>
#include <vector>

#include "io/file/netcdf_io.h"
#include "data_struct/spectra.h"
#include "data_struct/spectra_line.h"
#include "data_struct/scan_info.h"
#include "data_struct/params_override.h"

// ============================================================
// Layout constants — must stay in sync with netcdf_io.cpp
// ============================================================
static constexpr size_t NCDF_DIM2   = 1047808; // fixed read count per row
static constexpr size_t HDR_SZ      = 256;      // default header_size
static constexpr size_t MAX_DETS    = 4;        // MAX_NUM_SUPPORTED_DETECOTRS_PER_COL

// Sub-header offsets from the start of a column slot
static constexpr size_t RT_OFF = 32; // ELAPSED_REALTIME_OFFSET
static constexpr size_t LT_OFF = 34; // ELAPSED_LIVETIME_OFFSET
static constexpr size_t IC_OFF = 36; // INPUT_COUNTS_OFFSET
static constexpr size_t OC_OFF = 38; // OUTPUT_COUNTS_OFFSET

// ============================================================
// Encoding helpers
// ============================================================

// Convert seconds → two 16-bit halves (stored as floats) matching the
// hardware encoding read by _load_spectra.
static void encode_time(float secs, float& lo, float& hi)
{
    auto ii = static_cast<unsigned int>(secs / 320e-9f);
    lo = static_cast<float>(ii & 0x0000FFFFu);
    hi = static_cast<float>((ii >> 16) & 0x0000FFFFu);
}

static void encode_count(unsigned int total, float& lo, float& hi)
{
    lo = static_cast<float>(total & 0x0000FFFFu);
    hi = static_cast<float>((total >> 16) & 0x0000FFFFu);
}

// ============================================================
// NetCDF file creators
// ============================================================

struct SpectraFileParams
{
    size_t              detector    = 0;
    size_t              num_cols    = 2;
    size_t              spectra_sz  = 2048;
    float               livetime    = 0.10f;
    float               realtime    = 0.12f;
    float               input_rate  = 1000.0f; // counts/s
    float               output_rate = 900.0f;  // counts/s
    std::vector<float>  spectrum;               // spectra_sz elements
};

// Build a spectra NetCDF file in the exact binary layout that NetCDF_IO reads.
// Returns NC_NOERR on success.
static int make_spectra_nc(const std::string& path, const SpectraFileParams& p)
{
    const size_t idx_det  = p.detector % MAX_DETS;
    const size_t inc_sz   = HDR_SZ + p.spectra_sz * MAX_DETS;
    const size_t dim1_sz  = (p.detector > 3) ? 2u : 1u;
    const size_t dim1_idx = (p.detector > 3) ? 1u : 0u;

    // Flat buffer for one [dim1_sz][NCDF_DIM2] slice, zero-initialised.
    // We only write into the slice at dim1_idx.
    std::vector<float> buf(NCDF_DIM2, 0.0f);

    // --- main header ---
    buf[0]  = 21930.0f;
    buf[1]  = -21931.0f;
    buf[2]  = static_cast<float>(HDR_SZ);
    buf[8]  = static_cast<float>(p.num_cols);
    buf[12 + 2 * idx_det] = static_cast<float>(p.detector);
    buf[20] = static_cast<float>(p.spectra_sz);

    for (size_t m1 = 0; m1 < p.num_cols; ++m1)
    {
        size_t l = HDR_SZ + m1 * inc_sz;

        // sub-header magic
        buf[l]     = 13260.0f;
        buf[l + 1] = -13261.0f;

        float rt_lo, rt_hi, lt_lo, lt_hi, ic_lo, ic_hi, oc_lo, oc_hi;
        encode_time(p.realtime, rt_lo, rt_hi);
        encode_time(p.livetime, lt_lo, lt_hi);
        encode_count(static_cast<unsigned int>(p.input_rate  * p.livetime), ic_lo, ic_hi);
        encode_count(static_cast<unsigned int>(p.output_rate * p.realtime), oc_lo, oc_hi);

        // idx_det offsets the timing block by 8 elements per detector
        size_t det_off = idx_det * 8;
        buf[l + RT_OFF + det_off]     = rt_lo;
        buf[l + RT_OFF + det_off + 1] = rt_hi;
        buf[l + LT_OFF + det_off]     = lt_lo;
        buf[l + LT_OFF + det_off + 1] = lt_hi;
        buf[l + IC_OFF + det_off]     = ic_lo;
        buf[l + IC_OFF + det_off + 1] = ic_hi;
        buf[l + OC_OFF + det_off]     = oc_lo;
        buf[l + OC_OFF + det_off + 1] = oc_hi;

        // spectra data
        size_t spec_off = l + HDR_SZ + p.spectra_sz * idx_det;
        for (size_t k = 0; k < p.spectra_sz && k < p.spectrum.size(); ++k)
            buf[spec_off + k] = p.spectrum[k];
    }

    int ncid, varid, retval;
    int dim_ids[3];
    if ((retval = nc_create(path.c_str(), NC_CLOBBER, &ncid)) != NC_NOERR)
        return retval;

    if ((retval = nc_def_dim(ncid, "d0", 1,        &dim_ids[0])) != NC_NOERR) { nc_close(ncid); return retval; }
    if ((retval = nc_def_dim(ncid, "d1", dim1_sz,  &dim_ids[1])) != NC_NOERR) { nc_close(ncid); return retval; }
    if ((retval = nc_def_dim(ncid, "d2", NCDF_DIM2, &dim_ids[2])) != NC_NOERR) { nc_close(ncid); return retval; }
    if ((retval = nc_def_var(ncid, "array_data", NC_FLOAT, 3, dim_ids, &varid)) != NC_NOERR) { nc_close(ncid); return retval; }
    if ((retval = nc_enddef(ncid)) != NC_NOERR) { nc_close(ncid); return retval; }

    size_t start[] = {0, dim1_idx, 0};
    size_t count[] = {1, 1, NCDF_DIM2};
    retval = nc_put_vara_float(ncid, varid, start, count, buf.data());
    if (retval != NC_NOERR) { nc_close(ncid); return retval; }

    return nc_close(ncid);
}

// Build a scalers NetCDF file in the layout read by load_scalers_line.
// scaler_vals[scaler_idx][col_idx]  →  array_data[col_idx][0][scaler_idx]
static int make_scaler_nc(const std::string& path,
                           size_t num_cols, size_t num_scalers,
                           const std::vector<std::vector<float>>& scaler_vals)
{
    int ncid, varid, retval;
    int dim_ids[3];
    if ((retval = nc_create(path.c_str(), NC_CLOBBER, &ncid)) != NC_NOERR)
        return retval;

    if ((retval = nc_def_dim(ncid, "d0", num_cols,    &dim_ids[0])) != NC_NOERR) { nc_close(ncid); return retval; }
    if ((retval = nc_def_dim(ncid, "d1", 1,           &dim_ids[1])) != NC_NOERR) { nc_close(ncid); return retval; }
    if ((retval = nc_def_dim(ncid, "d2", num_scalers, &dim_ids[2])) != NC_NOERR) { nc_close(ncid); return retval; }
    if ((retval = nc_def_var(ncid, "array_data", NC_FLOAT, 3, dim_ids, &varid)) != NC_NOERR) { nc_close(ncid); return retval; }
    if ((retval = nc_enddef(ncid)) != NC_NOERR) { nc_close(ncid); return retval; }

    // Write each scaler: data[col][0][i] = scaler_vals[i][col]
    for (size_t i = 0; i < num_scalers; ++i)
    {
        size_t start[] = {0, 0, i};
        size_t count[] = {num_cols, 1, 1};
        retval = nc_put_vara_float(ncid, varid, start, count, scaler_vals[i].data());
        if (retval != NC_NOERR) { nc_close(ncid); return retval; }
    }
    return nc_close(ncid);
}

// ============================================================
// Test fixture
// ============================================================

class NetCDF_IO_Test : public ::testing::Test
{
protected:
    std::filesystem::path tmp_dir;

    void SetUp() override
    {
        tmp_dir = std::filesystem::temp_directory_path() / "xrf_netcdf_io_tests";
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

    static io::file::NetCDF_IO<float>* io()
    {
        return io::file::NetCDF_IO<float>::inst();
    }
};

// ============================================================
// load_spectra_line — error paths
// ============================================================

TEST_F(NetCDF_IO_Test, LoadSpectraLine_NonExistentFile_ReturnsZero)
{
    data_struct::Spectra_Line<float> line;
    line.resize_and_zero(2, 2048);
    EXPECT_EQ(io()->load_spectra_line("/no/such/path.nc", 0, &line), 0u);
}

TEST_F(NetCDF_IO_Test, LoadSpectraLine_InvalidMagicNumber_ReturnsZero)
{
    // File exists and has the right variable but all-zero content (magic missing).
    std::string path = tmp("bad_magic.nc");
    int ncid, varid;
    int dim_ids[3];
    ASSERT_EQ(NC_NOERR, nc_create(path.c_str(), NC_CLOBBER, &ncid));
    ASSERT_EQ(NC_NOERR, nc_def_dim(ncid, "d0", 1,        &dim_ids[0]));
    ASSERT_EQ(NC_NOERR, nc_def_dim(ncid, "d1", 1,        &dim_ids[1]));
    ASSERT_EQ(NC_NOERR, nc_def_dim(ncid, "d2", NCDF_DIM2, &dim_ids[2]));
    ASSERT_EQ(NC_NOERR, nc_def_var(ncid, "array_data", NC_FLOAT, 3, dim_ids, &varid));
    ASSERT_EQ(NC_NOERR, nc_enddef(ncid));
    nc_close(ncid);

    data_struct::Spectra_Line<float> line;
    line.resize_and_zero(2, 2048);
    EXPECT_EQ(io()->load_spectra_line(path, 0, &line), 0u);
}

TEST_F(NetCDF_IO_Test, LoadSpectraLine_WrongDetectorInHeader_ReturnsMaxSize)
{
    // File records detector 0; requesting detector 1 → header mismatch → (size_t)-1
    std::vector<float> spectrum(2048, 1.0f);
    SpectraFileParams p;
    p.detector = 0;  p.num_cols = 2;  p.spectra_sz = 2048;  p.spectrum = spectrum;
    p.livetime = 0.1f;  p.realtime = 0.1f;
    p.input_rate = 500.0f;  p.output_rate = 500.0f;

    std::string path = tmp("wrong_det.nc");
    ASSERT_EQ(NC_NOERR, make_spectra_nc(path, p));

    data_struct::Spectra_Line<float> line;
    line.resize_and_zero(2, 2048);
    EXPECT_EQ(io()->load_spectra_line(path, 1, &line), static_cast<size_t>(-1));
}

TEST_F(NetCDF_IO_Test, LoadSpectraLine_ExtendedDetectorRequiresDim1_2_ReturnsError)
{
    // Detectors 4-7 require dim1 == 2. A file built for detector 0 has dim1 == 1.
    std::vector<float> spectrum(2048, 0.0f);
    SpectraFileParams p;
    p.detector = 0;  p.num_cols = 2;  p.spectra_sz = 2048;  p.spectrum = spectrum;
    p.livetime = 0.1f;  p.realtime = 0.1f;
    p.input_rate = 500.0f;  p.output_rate = 500.0f;

    std::string path = tmp("det4_wrong_dim.nc");
    ASSERT_EQ(NC_NOERR, make_spectra_nc(path, p));

    data_struct::Spectra_Line<float> line;
    line.resize_and_zero(2, 2048);
    EXPECT_EQ(io()->load_spectra_line(path, 4, &line), static_cast<size_t>(-1));
}

// ============================================================
// load_spectra_line — success paths
// ============================================================

TEST_F(NetCDF_IO_Test, LoadSpectraLine_Detector0_ReturnsCorrectColumnCount)
{
    std::vector<float> spectrum(2048, 1.0f);
    SpectraFileParams p;
    p.detector = 0;  p.num_cols = 3;  p.spectra_sz = 2048;  p.spectrum = spectrum;
    p.livetime = 0.1f;  p.realtime = 0.1f;
    p.input_rate = 500.0f;  p.output_rate = 500.0f;

    std::string path = tmp("det0_count.nc");
    ASSERT_EQ(NC_NOERR, make_spectra_nc(path, p));

    data_struct::Spectra_Line<float> line;
    line.resize_and_zero(3, 2048);
    EXPECT_EQ(io()->load_spectra_line(path, 0, &line), 3u);
}

TEST_F(NetCDF_IO_Test, LoadSpectraLine_Detector0_LoadsCorrectSpectraValues)
{
    const size_t spectra_sz = 2048;
    const size_t num_cols   = 3;
    std::vector<float> spectrum(spectra_sz);
    for (size_t k = 0; k < spectra_sz; ++k)
        spectrum[k] = static_cast<float>(k + 1);

    SpectraFileParams p;
    p.detector = 0;  p.num_cols = num_cols;  p.spectra_sz = spectra_sz;  p.spectrum = spectrum;
    p.livetime = 0.1f;  p.realtime = 0.1f;
    p.input_rate = 500.0f;  p.output_rate = 500.0f;

    std::string path = tmp("det0_spectra.nc");
    ASSERT_EQ(NC_NOERR, make_spectra_nc(path, p));

    data_struct::Spectra_Line<float> line;
    line.resize_and_zero(num_cols, spectra_sz);
    ASSERT_EQ(io()->load_spectra_line(path, 0, &line), num_cols);

    for (size_t col = 0; col < num_cols; ++col)
        for (size_t k = 0; k < spectra_sz; ++k)
            EXPECT_FLOAT_EQ(line[col][k], spectrum[k])
                << "col=" << col << " k=" << k;
}

TEST_F(NetCDF_IO_Test, LoadSpectraLine_Detector0_HasCorrectRealtime)
{
    // Use equal input/output rates so recalc_elapsed_livetime gives realtime.
    const float rt = 0.12f;
    std::vector<float> spectrum(2048, 1.0f);
    SpectraFileParams p;
    p.detector = 0;  p.num_cols = 2;  p.spectra_sz = 2048;  p.spectrum = spectrum;
    p.livetime = rt;  p.realtime = rt;
    p.input_rate = 500.0f;  p.output_rate = 500.0f;

    std::string path = tmp("det0_timing.nc");
    ASSERT_EQ(NC_NOERR, make_spectra_nc(path, p));

    data_struct::Spectra_Line<float> line;
    line.resize_and_zero(2, 2048);
    ASSERT_EQ(io()->load_spectra_line(path, 0, &line), 2u);

    for (size_t col = 0; col < 2; ++col)
    {
        // Tolerance: one 320 ns tick
        EXPECT_NEAR(line[col].elapsed_realtime(), rt, 1e-4f) << "col=" << col;
        // recalc gives realtime * out_rate / in_rate = rt * 500/500 = rt
        EXPECT_NEAR(line[col].elapsed_livetime(), rt, 1e-4f) << "col=" << col;
    }
}

TEST_F(NetCDF_IO_Test, LoadSpectraLine_Detector1_LoadsCorrectSpectraValues)
{
    const size_t spectra_sz = 2048;
    std::vector<float> spectrum(spectra_sz, 42.0f);
    SpectraFileParams p;
    p.detector = 1;  p.num_cols = 2;  p.spectra_sz = spectra_sz;  p.spectrum = spectrum;
    p.livetime = 0.1f;  p.realtime = 0.1f;
    p.input_rate = 500.0f;  p.output_rate = 500.0f;

    std::string path = tmp("det1_spectra.nc");
    ASSERT_EQ(NC_NOERR, make_spectra_nc(path, p));

    data_struct::Spectra_Line<float> line;
    line.resize_and_zero(2, spectra_sz);
    ASSERT_EQ(io()->load_spectra_line(path, 1, &line), 2u);

    for (size_t col = 0; col < 2; ++col)
        for (size_t k = 0; k < spectra_sz; ++k)
            EXPECT_FLOAT_EQ(line[col][k], 42.0f) << "col=" << col << " k=" << k;
}

TEST_F(NetCDF_IO_Test, LoadSpectraLine_Detector4_LoadsCorrectSpectraValues)
{
    // Detectors 4-7 use the second dim1 slice.
    const size_t spectra_sz = 2048;
    std::vector<float> spectrum(spectra_sz);
    for (size_t k = 0; k < spectra_sz; ++k) spectrum[k] = static_cast<float>(k * 2);

    SpectraFileParams p;
    p.detector = 4;  p.num_cols = 2;  p.spectra_sz = spectra_sz;  p.spectrum = spectrum;
    p.livetime = 0.1f;  p.realtime = 0.1f;
    p.input_rate = 800.0f;  p.output_rate = 800.0f;

    std::string path = tmp("det4_spectra.nc");
    ASSERT_EQ(NC_NOERR, make_spectra_nc(path, p));

    data_struct::Spectra_Line<float> line;
    line.resize_and_zero(2, spectra_sz);
    ASSERT_EQ(io()->load_spectra_line(path, 4, &line), 2u);

    for (size_t col = 0; col < 2; ++col)
        for (size_t k = 0; k < spectra_sz; ++k)
            EXPECT_FLOAT_EQ(line[col][k], spectrum[k]) << "col=" << col << " k=" << k;
}

// ============================================================
// load_spectra_line_integrated
// ============================================================

TEST_F(NetCDF_IO_Test, LoadSpectraLineIntegrated_AccumulatesSpectraValues)
{
    const size_t spectra_sz = 2048;
    const size_t num_cols   = 2;
    std::vector<float> spectrum(spectra_sz, 10.0f);

    SpectraFileParams p;
    p.detector = 0;  p.num_cols = num_cols;  p.spectra_sz = spectra_sz;  p.spectrum = spectrum;
    p.livetime = 0.1f;  p.realtime = 0.1f;
    p.input_rate = 500.0f;  p.output_rate = 500.0f;

    std::string path = tmp("integrated.nc");
    ASSERT_EQ(NC_NOERR, make_spectra_nc(path, p));

    data_struct::Spectra<float> integrated(spectra_sz);
    integrated.setZero();
    ASSERT_EQ(io()->load_spectra_line_integrated(path, 0, num_cols, &integrated), num_cols);

    // Each of num_cols columns contributes 10.0 to every channel.
    for (size_t k = 0; k < spectra_sz; ++k)
        EXPECT_FLOAT_EQ(integrated(k), 10.0f * static_cast<float>(num_cols)) << "k=" << k;
}

TEST_F(NetCDF_IO_Test, LoadSpectraLineIntegrated_AccumulatesTiming)
{
    const size_t spectra_sz = 2048;
    const size_t num_cols   = 3;
    const float  livetime   = 0.10f;
    const float  realtime   = 0.12f;
    std::vector<float> spectrum(spectra_sz, 1.0f);

    SpectraFileParams p;
    p.detector = 0;  p.num_cols = num_cols;  p.spectra_sz = spectra_sz;  p.spectrum = spectrum;
    p.livetime = livetime;  p.realtime = realtime;
    p.input_rate = 1000.0f;  p.output_rate = 900.0f;

    std::string path = tmp("int_timing.nc");
    ASSERT_EQ(NC_NOERR, make_spectra_nc(path, p));

    data_struct::Spectra<float> integrated(spectra_sz);
    integrated.elapsed_livetime(0.0f);
    integrated.elapsed_realtime(0.0f);
    integrated.input_counts(0.0f);
    integrated.output_counts(0.0f);

    ASSERT_EQ(io()->load_spectra_line_integrated(path, 0, num_cols, &integrated), num_cols);

    // Timing is summed once per column, no recalc.
    EXPECT_NEAR(integrated.elapsed_livetime(), livetime * num_cols, 1e-4f);
    EXPECT_NEAR(integrated.elapsed_realtime(), realtime * num_cols, 1e-4f);
}

// ============================================================
// load_spectra_line_with_callback
// ============================================================

TEST_F(NetCDF_IO_Test, LoadSpectraLineWithCallback_InvalidFile_ReturnsFalse)
{
    data_struct::IO_Callback_Func_Def<float> cb =
        [](size_t, size_t, size_t, size_t, size_t, data_struct::Spectra<float>*, void*) {};

    EXPECT_FALSE(io()->load_spectra_line_with_callback(
        "/no/such/file.nc", {0}, 0, 1, 10, cb, nullptr));
}

TEST_F(NetCDF_IO_Test, LoadSpectraLineWithCallback_InvokedOncePerColumn)
{
    const size_t spectra_sz = 2048;
    const size_t num_cols   = 3;
    std::vector<float> spectrum(spectra_sz);
    for (size_t k = 0; k < spectra_sz; ++k) spectrum[k] = static_cast<float>(k + 1);

    SpectraFileParams p;
    p.detector = 0;  p.num_cols = num_cols;  p.spectra_sz = spectra_sz;  p.spectrum = spectrum;
    p.livetime = 0.1f;  p.realtime = 0.1f;
    p.input_rate = 500.0f;  p.output_rate = 500.0f;

    std::string path = tmp("callback.nc");
    ASSERT_EQ(NC_NOERR, make_spectra_nc(path, p));

    size_t call_count = 0;
    std::vector<size_t> col_indices;
    std::vector<float>  first_channels;

    data_struct::IO_Callback_Func_Def<float> cb =
        [&](size_t /*row*/, size_t col, size_t /*max_rows*/, size_t /*max_cols*/,
            size_t /*det*/, data_struct::Spectra<float>* spec, void*)
        {
            ++call_count;
            col_indices.push_back(col);
            if (spec) first_channels.push_back((*spec)(0));
        };

    EXPECT_TRUE(io()->load_spectra_line_with_callback(path, {0}, 0, 1, num_cols, cb, nullptr));
    ASSERT_EQ(call_count, num_cols);

    for (size_t i = 0; i < num_cols; ++i)
        EXPECT_EQ(col_indices[i], i) << "i=" << i;

    for (float ch : first_channels)
        EXPECT_FLOAT_EQ(ch, 1.0f); // first bin is 1
}

TEST_F(NetCDF_IO_Test, LoadSpectraLineWithCallback_CorrectDetectorIdPassedToCallback)
{
    std::vector<float> spectrum(2048, 5.0f);
    SpectraFileParams p;
    p.detector = 0;  p.num_cols = 2;  p.spectra_sz = 2048;  p.spectrum = spectrum;
    p.livetime = 0.1f;  p.realtime = 0.1f;
    p.input_rate = 500.0f;  p.output_rate = 500.0f;

    std::string path = tmp("cb_det.nc");
    ASSERT_EQ(NC_NOERR, make_spectra_nc(path, p));

    std::vector<size_t> det_ids;
    data_struct::IO_Callback_Func_Def<float> cb =
        [&](size_t, size_t, size_t, size_t, size_t det, data_struct::Spectra<float>*, void*)
        { det_ids.push_back(det); };

    EXPECT_TRUE(io()->load_spectra_line_with_callback(path, {0}, 0, 1, 2, cb, nullptr));
    ASSERT_EQ(det_ids.size(), 2u);
    for (size_t id : det_ids) EXPECT_EQ(id, 0u);
}

TEST_F(NetCDF_IO_Test, LoadSpectraLineWithCallback_CorrectSpectraDataPassedToCallback)
{
    const size_t spectra_sz = 2048;
    std::vector<float> spectrum(spectra_sz);
    for (size_t k = 0; k < spectra_sz; ++k) spectrum[k] = static_cast<float>(k + 1);

    SpectraFileParams p;
    p.detector = 0;  p.num_cols = 2;  p.spectra_sz = spectra_sz;  p.spectrum = spectrum;
    p.livetime = 0.1f;  p.realtime = 0.1f;
    p.input_rate = 500.0f;  p.output_rate = 500.0f;

    std::string path = tmp("cb_spec.nc");
    ASSERT_EQ(NC_NOERR, make_spectra_nc(path, p));

    std::vector<data_struct::Spectra<float>> received;
    data_struct::IO_Callback_Func_Def<float> cb =
        [&](size_t, size_t, size_t, size_t, size_t, data_struct::Spectra<float>* spec, void*)
        { if (spec) received.push_back(*spec); };

    EXPECT_TRUE(io()->load_spectra_line_with_callback(path, {0}, 0, 1, 2, cb, nullptr));
    ASSERT_EQ(received.size(), 2u);
    for (auto& spec : received)
        for (size_t k = 0; k < spectra_sz; ++k)
            EXPECT_FLOAT_EQ(spec(k), spectrum[k]) << "k=" << k;
}

// ============================================================
// load_scalers_line — error paths
// ============================================================

TEST_F(NetCDF_IO_Test, LoadScalersLine_NullScanInfo_ReturnsZero)
{
    EXPECT_EQ(io()->load_scalers_line("any.nc", "s", 0, nullptr, nullptr), 0u);
}

TEST_F(NetCDF_IO_Test, LoadScalersLine_EmptyScalerMaps_ReturnsZero)
{
    data_struct::Scan_Info<float> info; // scaler_maps empty
    EXPECT_EQ(io()->load_scalers_line("any.nc", "s", 0, &info, nullptr), 0u);
}

TEST_F(NetCDF_IO_Test, LoadScalersLine_RowOutOfBounds_ReturnsZero)
{
    data_struct::Scan_Info<float> info;
    info.initialize_scaler_map_with_dims("sc", 1 /*rows*/, 2 /*cols*/, false, "s0");
    // Row 5 requested but only 1 row allocated.
    EXPECT_EQ(io()->load_scalers_line("any.nc", "s", 5, &info, nullptr), 0u);
}

// ============================================================
// load_scalers_line — success paths
// ============================================================

TEST_F(NetCDF_IO_Test, LoadScalersLine_ValidFile_LoadsCorrectValues)
{
    const size_t num_cols    = 4;
    const size_t num_scalers = 2;
    const size_t row         = 0;

    // scaler_vals[scaler_idx][col_idx]
    std::vector<std::vector<float>> vals = {
        {1.0f, 2.0f, 3.0f, 4.0f}, // scaler index 0
        {5.0f, 6.0f, 7.0f, 8.0f}, // scaler index 1
    };

    std::string path = tmp("scalers.nc");
    ASSERT_EQ(NC_NOERR, make_scaler_nc(path, num_cols, num_scalers, vals));

    data_struct::Scan_Info<float> info;
    // The scaler map unit must equal tag + std::to_string(scaler_index).
    info.initialize_scaler_map_with_dims("scaler_a", 1, num_cols, false, "myscaler0");
    info.initialize_scaler_map_with_dims("scaler_b", 1, num_cols, false, "myscaler1");

    size_t loaded = io()->load_scalers_line(path, "myscaler", row, &info, nullptr);
    EXPECT_EQ(loaded, num_cols);

    auto& a = info.scaler_maps.at("scaler_a");
    auto& b = info.scaler_maps.at("scaler_b");
    for (size_t j = 0; j < num_cols; ++j)
    {
        EXPECT_FLOAT_EQ(a->values(row, j), vals[0][j]) << "j=" << j;
        EXPECT_FLOAT_EQ(b->values(row, j), vals[1][j]) << "j=" << j;
    }
}

TEST_F(NetCDF_IO_Test, LoadScalersLine_WithScalingFactor_AppliesMultiplier)
{
    const size_t num_cols    = 3;
    const size_t num_scalers = 1;
    const size_t row         = 0;
    const float  factor      = 2.5f;

    std::vector<std::vector<float>> vals = {{10.0f, 20.0f, 30.0f}};

    std::string path = tmp("scalers_factor.nc");
    ASSERT_EQ(NC_NOERR, make_scaler_nc(path, num_cols, num_scalers, vals));

    data_struct::Scan_Info<float> info;
    info.initialize_scaler_map_with_dims("my_scaler", 1, num_cols, false, "tag0");

    data_struct::Params_Override<float> params;
    params.scaling_factors["my_scaler"] = factor;

    size_t loaded = io()->load_scalers_line(path, "tag", row, &info, &params);
    EXPECT_EQ(loaded, num_cols);

    auto& m = info.scaler_maps.at("my_scaler");
    for (size_t j = 0; j < num_cols; ++j)
        EXPECT_FLOAT_EQ(m->values(row, j), vals[0][j] * factor) << "j=" << j;
}

TEST_F(NetCDF_IO_Test, LoadScalersLine_NonZeroRow_LoadsIntoCorrectRow)
{
    const size_t num_cols    = 2;
    const size_t num_scalers = 1;

    std::vector<std::vector<float>> vals = {{7.0f, 8.0f}};

    std::string path = tmp("scalers_row.nc");
    ASSERT_EQ(NC_NOERR, make_scaler_nc(path, num_cols, num_scalers, vals));

    data_struct::Scan_Info<float> info;
    // Allocate 3 rows so row=2 is valid.
    info.initialize_scaler_map_with_dims("sc", 3, num_cols, false, "pfx0");

    size_t loaded = io()->load_scalers_line(path, "pfx", 2 /*row*/, &info, nullptr);
    EXPECT_EQ(loaded, num_cols);

    auto& m = info.scaler_maps.at("sc");
    // Rows 0 and 1 should remain zero; row 2 should have the loaded values.
    EXPECT_FLOAT_EQ(m->values(0, 0), 0.0f);
    EXPECT_FLOAT_EQ(m->values(2, 0), vals[0][0]);
    EXPECT_FLOAT_EQ(m->values(2, 1), vals[0][1]);
}
