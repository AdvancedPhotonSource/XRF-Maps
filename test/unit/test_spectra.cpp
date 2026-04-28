// Unit tests for data_struct::Spectra<_T> and related free functions in spectra.h
//
// Covers:
//   - Spectra constructors (default, size, size+timing, copy, ArrayTr, ArrayTr+timing)
//   - Accessors (elapsed_livetime, elapsed_realtime, input_counts, output_counts)
//   - recalc_elapsed_livetime()
//   - set_nan_to_near_zero()
//   - add(const Spectra<_T>&)        – same-type, by value
//   - add(const Spectra<T_real>*)    – cross-type, by pointer
//   - sub_spectra()
//   - convolve1d() – both overloads
//   - snip_background()
//
// Build with BUILD_TESTS=ON:
//   cmake -DBUILD_TESTS=ON ..
//   cmake --build . --target test_spectra
//   ctest  (or ./bin/test_spectra)

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include "data_struct/spectra.h"

using namespace data_struct;

// ============================================================
// Constructors
// ============================================================

TEST(SpectraConstructors, DefaultConstructorFloat)
{
    Spectra<float> s;
    EXPECT_FLOAT_EQ(s.elapsed_livetime(), (float)default_time_and_io_counts);
    EXPECT_FLOAT_EQ(s.elapsed_realtime(), (float)default_time_and_io_counts);
    EXPECT_FLOAT_EQ(s.input_counts(),  0.0f);
    EXPECT_FLOAT_EQ(s.output_counts(), 0.0f);
}

TEST(SpectraConstructors, DefaultConstructorDouble)
{
    Spectra<double> s;
    EXPECT_DOUBLE_EQ(s.elapsed_livetime(), (double)default_time_and_io_counts);
    EXPECT_DOUBLE_EQ(s.elapsed_realtime(), (double)default_time_and_io_counts);
    EXPECT_DOUBLE_EQ(s.input_counts(),  0.0);
    EXPECT_DOUBLE_EQ(s.output_counts(), 0.0);
}

TEST(SpectraConstructors, SizeConstructorZeroInitialised)
{
    Spectra<float> s(8);
    EXPECT_EQ(s.size(), 8);
    for (Eigen::Index i = 0; i < s.size(); ++i)
        EXPECT_FLOAT_EQ(s[i], 0.0f);
    EXPECT_FLOAT_EQ(s.elapsed_livetime(), (float)default_time_and_io_counts);
    EXPECT_FLOAT_EQ(s.elapsed_realtime(), (float)default_time_and_io_counts);
    EXPECT_FLOAT_EQ(s.input_counts(),  0.0f);
    EXPECT_FLOAT_EQ(s.output_counts(), 0.0f);
}

TEST(SpectraConstructors, SizeWithTimingConstructor)
{
    Spectra<double> s(6, 1.5, 2.0, 100.0, 95.0);
    EXPECT_EQ(s.size(), 6);
    for (Eigen::Index i = 0; i < s.size(); ++i)
        EXPECT_DOUBLE_EQ(s[i], 0.0);
    EXPECT_DOUBLE_EQ(s.elapsed_livetime(), 1.5);
    EXPECT_DOUBLE_EQ(s.elapsed_realtime(), 2.0);
    EXPECT_DOUBLE_EQ(s.input_counts(),    100.0);
    EXPECT_DOUBLE_EQ(s.output_counts(),    95.0);
}

TEST(SpectraConstructors, CopyConstructorPreservesDataAndTiming)
{
    Spectra<float> orig(4, 1.0f, 2.0f, 50.0f, 45.0f);
    orig[0] = 1.0f; orig[1] = 2.0f; orig[2] = 3.0f; orig[3] = 4.0f;

    Spectra<float> copy(orig);

    EXPECT_EQ(copy.size(), 4);
    for (Eigen::Index i = 0; i < 4; ++i)
        EXPECT_FLOAT_EQ(copy[i], orig[i]);
    EXPECT_FLOAT_EQ(copy.elapsed_livetime(), 1.0f);
    EXPECT_FLOAT_EQ(copy.elapsed_realtime(), 2.0f);
    EXPECT_FLOAT_EQ(copy.input_counts(),    50.0f);
    EXPECT_FLOAT_EQ(copy.output_counts(),   45.0f);
}

TEST(SpectraConstructors, CopyIsIndependent)
{
    Spectra<float> orig(3);
    orig[0] = 1.0f; orig[1] = 2.0f; orig[2] = 3.0f;
    Spectra<float> copy(orig);

    copy[0] = 99.0f;
    EXPECT_FLOAT_EQ(orig[0], 1.0f);  // orig unchanged
}

TEST(SpectraConstructors, FromArrayDefaultTiming)
{
    ArrayTr<float> arr(5);
    arr << 10.0f, 20.0f, 30.0f, 40.0f, 50.0f;

    Spectra<float> s(arr);

    EXPECT_EQ(s.size(), 5);
    for (Eigen::Index i = 0; i < s.size(); ++i)
        EXPECT_FLOAT_EQ(s[i], arr[i]);
    EXPECT_FLOAT_EQ(s.elapsed_livetime(), (float)default_time_and_io_counts);
    EXPECT_FLOAT_EQ(s.elapsed_realtime(), (float)default_time_and_io_counts);
}

TEST(SpectraConstructors, FromArrayWithTiming)
{
    ArrayTr<double> arr(3);
    arr << 10.0, 20.0, 30.0;

    Spectra<double> s(arr, 0.5, 1.0, 200.0, 180.0);

    EXPECT_EQ(s.size(), 3);
    EXPECT_DOUBLE_EQ(s[0], 10.0);
    EXPECT_DOUBLE_EQ(s[1], 20.0);
    EXPECT_DOUBLE_EQ(s[2], 30.0);
    EXPECT_DOUBLE_EQ(s.elapsed_livetime(),  0.5);
    EXPECT_DOUBLE_EQ(s.elapsed_realtime(),  1.0);
    EXPECT_DOUBLE_EQ(s.input_counts(),    200.0);
    EXPECT_DOUBLE_EQ(s.output_counts(),   180.0);
}

// ============================================================
// Accessors
// ============================================================

TEST(SpectraAccessors, SetAndGet)
{
    Spectra<float> s(2);
    s.elapsed_livetime(1.1f);
    s.elapsed_realtime(2.2f);
    s.input_counts(300.0f);
    s.output_counts(290.0f);

    EXPECT_FLOAT_EQ(s.elapsed_livetime(),  1.1f);
    EXPECT_FLOAT_EQ(s.elapsed_realtime(),  2.2f);
    EXPECT_FLOAT_EQ(s.input_counts(),    300.0f);
    EXPECT_FLOAT_EQ(s.output_counts(),   290.0f);
}

// ============================================================
// recalc_elapsed_livetime
// ============================================================

TEST(RecalcLivetime, NonZeroCounts)
{
    Spectra<double> s(2, 0.0, 2.0, 1000.0, 500.0);
    s.recalc_elapsed_livetime();
    // livetime = realtime * output / input = 2.0 * 500.0 / 1000.0 = 1.0
    EXPECT_DOUBLE_EQ(s.elapsed_livetime(), 1.0);
}

TEST(RecalcLivetime, ZeroInputCountsFallsBackToRealtime)
{
    Spectra<double> s(2, 0.0, 3.0, 0.0, 500.0);
    s.recalc_elapsed_livetime();
    EXPECT_DOUBLE_EQ(s.elapsed_livetime(), 3.0);
}

TEST(RecalcLivetime, ZeroOutputCountsFallsBackToRealtime)
{
    Spectra<float> s(2, 0.0f, 3.0f, 500.0f, 0.0f);
    s.recalc_elapsed_livetime();
    EXPECT_FLOAT_EQ(s.elapsed_livetime(), 3.0f);
}

TEST(RecalcLivetime, BothCountsZeroFallsBackToRealtime)
{
    Spectra<float> s(2, 0.0f, 5.0f, 0.0f, 0.0f);
    s.recalc_elapsed_livetime();
    EXPECT_FLOAT_EQ(s.elapsed_livetime(), 5.0f);
}

TEST(RecalcLivetime, EqualInputOutputGivesRealtime)
{
    Spectra<double> s(2, 0.0, 4.0, 200.0, 200.0);
    s.recalc_elapsed_livetime();
    EXPECT_DOUBLE_EQ(s.elapsed_livetime(), 4.0);
}

// ============================================================
// set_nan_to_near_zero
// ============================================================

TEST(SetNanToNearZero, NanReplaced)
{
    Spectra<float> s(4);
    s[0] = 1.0f;
    s[1] = std::numeric_limits<float>::quiet_NaN();
    s[2] = 3.0f;
    s[3] = 4.0f;

    s.set_nan_to_near_zero();

    EXPECT_FLOAT_EQ(s[0], 1.0f);
    EXPECT_FLOAT_EQ(s[1], (float)default_time_and_io_counts);
    EXPECT_FLOAT_EQ(s[2], 3.0f);
    EXPECT_FLOAT_EQ(s[3], 4.0f);
}

TEST(SetNanToNearZero, PositiveInfReplaced)
{
    Spectra<double> s(2);
    s[0] = std::numeric_limits<double>::infinity();
    s[1] = 5.0;

    s.set_nan_to_near_zero();

    EXPECT_DOUBLE_EQ(s[0], (double)default_time_and_io_counts);
    EXPECT_DOUBLE_EQ(s[1], 5.0);
}

TEST(SetNanToNearZero, NegativeInfReplaced)
{
    Spectra<float> s(2);
    s[0] = -std::numeric_limits<float>::infinity();
    s[1] = 2.0f;

    s.set_nan_to_near_zero();

    EXPECT_FLOAT_EQ(s[0], (float)default_time_and_io_counts);
    EXPECT_FLOAT_EQ(s[1], 2.0f);
}

TEST(SetNanToNearZero, AllFiniteUnchanged)
{
    Spectra<float> s(3);
    s[0] = 1.0f; s[1] = 2.0f; s[2] = 3.0f;

    s.set_nan_to_near_zero();

    EXPECT_FLOAT_EQ(s[0], 1.0f);
    EXPECT_FLOAT_EQ(s[1], 2.0f);
    EXPECT_FLOAT_EQ(s[2], 3.0f);
}

// ============================================================
// add – same type, by value
// ============================================================

TEST(SpectraAddSameType, ArrayValuesAndTimingAccumulate)
{
    Spectra<float> a(3, 1.0f, 2.0f, 100.0f, 90.0f);
    a[0] = 1.0f; a[1] = 2.0f; a[2] = 3.0f;

    Spectra<float> b(3, 0.5f, 1.0f, 50.0f, 45.0f);
    b[0] = 4.0f; b[1] = 5.0f; b[2] = 6.0f;

    EXPECT_TRUE(a.add(b));

    EXPECT_FLOAT_EQ(a[0], 5.0f);
    EXPECT_FLOAT_EQ(a[1], 7.0f);
    EXPECT_FLOAT_EQ(a[2], 9.0f);
    EXPECT_FLOAT_EQ(a.elapsed_livetime(),  1.5f);
    EXPECT_FLOAT_EQ(a.elapsed_realtime(),  3.0f);
    EXPECT_FLOAT_EQ(a.input_counts(),    150.0f);
    EXPECT_FLOAT_EQ(a.output_counts(),   135.0f);
}

TEST(SpectraAddSameType, SizeMismatchReturnsFalse)
{
    Spectra<double> a(3);
    Spectra<double> b(5);
    EXPECT_FALSE(a.add(b));
}

TEST(SpectraAddSameType, SizeMismatchDoesNotModify)
{
    Spectra<float> a(3, 1.0f, 2.0f, 100.0f, 90.0f);
    a[0] = 1.0f; a[1] = 2.0f; a[2] = 3.0f;
    Spectra<float> b(5);

    a.add(b);

    EXPECT_FLOAT_EQ(a[0], 1.0f);
    EXPECT_FLOAT_EQ(a.elapsed_livetime(), 1.0f);
}

TEST(SpectraAddSameType, NonFiniteLivetimeSkipped)
{
    Spectra<float> a(2, 1.0f, 2.0f, 100.0f, 90.0f);
    a[0] = 1.0f; a[1] = 1.0f;

    Spectra<float> b(2, std::numeric_limits<float>::quiet_NaN(), 0.5f, 10.0f, 9.0f);
    b[0] = 0.0f; b[1] = 0.0f;

    a.add(b);

    EXPECT_FLOAT_EQ(a.elapsed_livetime(), 1.0f);  // NaN was skipped
    EXPECT_FLOAT_EQ(a.elapsed_realtime(), 2.5f);
}

TEST(SpectraAddSameType, NonFiniteRealtimeSkipped)
{
    Spectra<double> a(2, 1.0, 2.0, 100.0, 90.0);
    a[0] = 1.0; a[1] = 1.0;

    Spectra<double> b(2, 0.5, std::numeric_limits<double>::infinity(), 10.0, 9.0);
    b[0] = 0.0; b[1] = 0.0;

    a.add(b);

    EXPECT_DOUBLE_EQ(a.elapsed_realtime(), 2.0);  // inf was skipped
    EXPECT_DOUBLE_EQ(a.elapsed_livetime(), 1.5);
}

TEST(SpectraAddSameType, NonFiniteInputCountsSkipped)
{
    Spectra<float> a(2, 1.0f, 2.0f, 100.0f, 90.0f);
    a[0] = 1.0f; a[1] = 1.0f;

    Spectra<float> b(2, 0.5f, 0.5f, std::numeric_limits<float>::quiet_NaN(), 9.0f);
    b[0] = 0.0f; b[1] = 0.0f;

    a.add(b);

    EXPECT_FLOAT_EQ(a.input_counts(), 100.0f);  // NaN was skipped
    EXPECT_FLOAT_EQ(a.output_counts(), 99.0f);
}

TEST(SpectraAddSameType, NonFiniteOutputCountsSkipped)
{
    Spectra<float> a(2, 1.0f, 2.0f, 100.0f, 90.0f);
    a[0] = 1.0f; a[1] = 1.0f;

    Spectra<float> b(2, 0.5f, 0.5f, 10.0f, std::numeric_limits<float>::infinity());
    b[0] = 0.0f; b[1] = 0.0f;

    a.add(b);

    EXPECT_FLOAT_EQ(a.output_counts(), 90.0f);  // inf was skipped
    EXPECT_FLOAT_EQ(a.input_counts(), 110.0f);
}

TEST(SpectraAddSameType, DoubleAccumulation)
{
    Spectra<double> acc(2, 0.0, 0.0, 0.0, 0.0);
    acc[0] = 0.0; acc[1] = 0.0;

    for (int i = 0; i < 5; ++i)
    {
        Spectra<double> s(2, 1.0, 2.0, 100.0, 90.0);
        s[0] = 10.0; s[1] = 20.0;
        acc.add(s);
    }

    EXPECT_DOUBLE_EQ(acc[0], 50.0);
    EXPECT_DOUBLE_EQ(acc[1], 100.0);
    EXPECT_DOUBLE_EQ(acc.elapsed_livetime(), 5.0);
    EXPECT_DOUBLE_EQ(acc.elapsed_realtime(), 10.0);
    EXPECT_DOUBLE_EQ(acc.input_counts(),  500.0);
    EXPECT_DOUBLE_EQ(acc.output_counts(), 450.0);
}

// ============================================================
// add – cross-type, by pointer
// ============================================================

TEST(SpectraAddCrossType, NullptrReturnsFalse)
{
    Spectra<float> a(4);
    const Spectra<double>* p = nullptr;
    EXPECT_FALSE(a.add(p));
}

TEST(SpectraAddCrossType, SizeMismatchReturnsFalse)
{
    Spectra<float> a(3);
    Spectra<double> b(5);
    EXPECT_FALSE(a.add(&b));
}

TEST(SpectraAddCrossType, TimingAccumulatedFromZero)
{
    // Start with a = 0 so any doubling artefact in array values is hidden.
    Spectra<float> a(2, 1.0f, 2.0f, 100.0f, 90.0f);
    a[0] = 0.0f; a[1] = 0.0f;

    Spectra<double> b(2, 0.5, 1.0, 50.0, 45.0);
    b[0] = 0.0; b[1] = 0.0;

    EXPECT_TRUE(a.add(&b));

    EXPECT_FLOAT_EQ(a.elapsed_livetime(),  1.5f);
    EXPECT_FLOAT_EQ(a.elapsed_realtime(),  3.0f);
    EXPECT_FLOAT_EQ(a.input_counts(),    150.0f);
    EXPECT_FLOAT_EQ(a.output_counts(),   135.0f);
}

TEST(SpectraAddCrossType, NonFiniteLivetimeSkipped)
{
    Spectra<float> a(2, 1.0f, 2.0f, 100.0f, 90.0f);
    a[0] = 0.0f; a[1] = 0.0f;

    Spectra<double> b(2, std::numeric_limits<double>::quiet_NaN(), 0.5, 10.0, 9.0);
    b[0] = 0.0; b[1] = 0.0;

    a.add(&b);

    EXPECT_FLOAT_EQ(a.elapsed_livetime(), 1.0f);  // NaN was skipped
    EXPECT_FLOAT_EQ(a.elapsed_realtime(), 2.5f);
}

// ============================================================
// sub_spectra
// ============================================================

TEST(SubSpectra, MiddleSegment)
{
    Spectra<float> s(6, 1.0f, 2.0f, 100.0f, 90.0f);
    s[0] = 10.0f; s[1] = 20.0f; s[2] = 30.0f;
    s[3] = 40.0f; s[4] = 50.0f; s[5] = 60.0f;

    Spectra<float> sub = s.sub_spectra(2, 3);

    EXPECT_EQ(sub.size(), 3);
    EXPECT_FLOAT_EQ(sub[0], 30.0f);
    EXPECT_FLOAT_EQ(sub[1], 40.0f);
    EXPECT_FLOAT_EQ(sub[2], 50.0f);
}

TEST(SubSpectra, InheritsTiming)
{
    Spectra<double> s(5, 1.5, 3.0, 200.0, 180.0);
    for (Eigen::Index i = 0; i < 5; ++i)
        s[i] = (double)i;

    Spectra<double> sub = s.sub_spectra(1, 3);

    EXPECT_DOUBLE_EQ(sub.elapsed_livetime(),  1.5);
    EXPECT_DOUBLE_EQ(sub.elapsed_realtime(),  3.0);
    EXPECT_DOUBLE_EQ(sub.input_counts(),    200.0);
    EXPECT_DOUBLE_EQ(sub.output_counts(),   180.0);
}

TEST(SubSpectra, FullRangeMatches)
{
    Spectra<float> s(4);
    s[0] = 1.0f; s[1] = 2.0f; s[2] = 3.0f; s[3] = 4.0f;

    Spectra<float> sub = s.sub_spectra(0, 4);

    EXPECT_EQ(sub.size(), 4);
    for (Eigen::Index i = 0; i < 4; ++i)
        EXPECT_FLOAT_EQ(sub[i], s[i]);
}

TEST(SubSpectra, SingleElement)
{
    Spectra<double> s(5);
    for (Eigen::Index i = 0; i < 5; ++i)
        s[i] = (double)(i * 10);

    Spectra<double> sub = s.sub_spectra(3, 1);

    EXPECT_EQ(sub.size(), 1);
    EXPECT_DOUBLE_EQ(sub[0], 30.0);
}

// ============================================================
// convolve1d
// ============================================================

TEST(Convolve1d, ConstantSignalPreservedInInterior)
{
    // A constant array convolved with a uniform boxcar normalised by its size
    // should return the same constant in the valid interior positions.
    ArrayTr<float> arr(10);
    arr.setConstant(9.0f);
    ArrayTr<float> boxcar(3);
    boxcar.setConstant(1.0f);

    ArrayTr<float> result = convolve1d(arr, boxcar);

    EXPECT_EQ(result.size(), arr.size());
    for (Eigen::Index i = 1; i < result.size() - 1; ++i)
        EXPECT_NEAR(result[i], 9.0f, 1e-5f);
}

TEST(Convolve1d, SizeOverloadMatchesExplicit)
{
    ArrayTr<double> arr(8);
    arr.setConstant(4.0);

    ArrayTr<double> boxcar(3);
    boxcar.setConstant(1.0);

    ArrayTr<double> r1 = convolve1d(arr, boxcar);
    ArrayTr<double> r2 = convolve1d(arr, (size_t)3);

    ASSERT_EQ(r1.size(), r2.size());
    for (Eigen::Index i = 0; i < r1.size(); ++i)
        EXPECT_NEAR(r1[i], r2[i], 1e-10);
}

TEST(Convolve1d, OutputSizeMatchesInput)
{
    ArrayTr<float> arr(20);
    arr.setConstant(1.0f);

    ArrayTr<float> result = convolve1d(arr, (size_t)5);

    EXPECT_EQ(result.size(), (Eigen::Index)20);
}

// ============================================================
// snip_background
// ============================================================

TEST(SnipBackground, NullptrReturnsEmpty)
{
    ArrayTr<float> bg = snip_background<float>(nullptr, 0.0f, 0.01f, 0.0f, 0.5f, 0.0f, 63.0f);
    EXPECT_EQ(bg.size(), 0);
}

TEST(SnipBackground, EmptySpectraReturnsEmpty)
{
    Spectra<float> s(0);
    ArrayTr<float> bg = snip_background<float>(&s, 0.0f, 0.01f, 0.0f, 0.5f, 0.0f, 0.0f);
    EXPECT_EQ(bg.size(), 0);
}

TEST(SnipBackground, OutputSizeMatchesInput)
{
    Spectra<float> s(64);
    for (Eigen::Index i = 0; i < s.size(); ++i)
        s[i] = (float)(i + 1);

    ArrayTr<float> bg = snip_background<float>(&s, 0.0f, 0.01f, 0.0f, 1.0f, 0.0f, 63.0f);

    EXPECT_EQ(bg.size(), s.size());
}

TEST(SnipBackground, BackgroundNeverExceedsSignal)
{
    // SNIP is a peak-stripping algorithm: the estimated background must be
    // <= the original signal at every channel.
    Spectra<float> s(64);
    s.setConstant(10.0f);
    s[30] = 1000.0f;  // sharp peak

    ArrayTr<float> bg = snip_background<float>(&s, 0.0f, 0.01f, 0.0f, 1.0f, 0.0f, 63.0f);

    for (Eigen::Index i = 0; i < bg.size(); ++i)
        EXPECT_LE(bg[i], s[i] + 1e-3f) << "background exceeds signal at channel " << i;
}

TEST(SnipBackground, OutputIsAllFinite)
{
    Spectra<double> s(32);
    for (Eigen::Index i = 0; i < s.size(); ++i)
        s[i] = (double)(i + 1) * 10.0;

    ArrayTr<double> bg = snip_background<double>(&s, 0.0, 0.01, 0.0, 1.0, 0.0, 31.0);

    for (Eigen::Index i = 0; i < bg.size(); ++i)
        EXPECT_TRUE(std::isfinite(bg[i])) << "bg[" << i << "] is not finite";
}

TEST(SnipBackground, BackgroundNonNegative)
{
    Spectra<float> s(32);
    for (Eigen::Index i = 0; i < s.size(); ++i)
        s[i] = (float)(i + 1);

    ArrayTr<float> bg = snip_background<float>(&s, 0.0f, 0.01f, 0.0f, 1.0f, 0.0f, 31.0f);

    for (Eigen::Index i = 0; i < bg.size(); ++i)
        EXPECT_GE(bg[i], 0.0f) << "background is negative at channel " << i;
}
