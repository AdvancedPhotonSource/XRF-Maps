#include <gtest/gtest.h>

#include <cstdint>
#include <cstring>
#include <memory>
#include <string>

#include "io/net/basic_serializer.h"

namespace
{

using Serializer = io::net::Basic_Serializer<float>;
using StreamBlock = data_struct::Stream_Block<float>;
using StreamBlockPtr = data_struct::Stream_Block_Ptr<float>;

StreamBlockPtr make_block()
{
    auto block = std::make_shared<StreamBlock>(2, 3, 4, 5, 6, "dataset/dir", "scan.mda");
    block->theta = 12.5f;
    block->spectra = new data_struct::Spectra<float>(4);
    block->spectra->setZero();
    (*block->spectra)[1] = 10.0f;
    (*block->spectra)[3] = 20.0f;
    block->fitting_blocks[data_struct::Fitting_Routines::ROI].fit_counts["Fe"] = 42.0f;
    return block;
}

TEST(BasicSerializer, SpectraRoundTrip)
{
    Serializer serializer;
    const std::string encoded = serializer.encode_spectra(make_block());
    auto decoded = serializer.decode_spectra(const_cast<char*>(encoded.data()), encoded.size());

    ASSERT_NE(decoded, nullptr);
    ASSERT_NE(decoded->spectra, nullptr);
    EXPECT_EQ(decoded->detector_number(), 2);
    EXPECT_EQ(*decoded->dataset_name, "scan.mda");
    EXPECT_EQ(*decoded->dataset_directory, "dataset/dir");
    ASSERT_EQ(decoded->spectra->size(), 4);
    EXPECT_FLOAT_EQ((*decoded->spectra)[1], 10.0f);
    EXPECT_FLOAT_EQ((*decoded->spectra)[3], 20.0f);
}

TEST(BasicSerializer, RejectsEveryTruncatedSpectraPacket)
{
    Serializer serializer;
    const std::string encoded = serializer.encode_spectra(make_block());

    for (size_t length = 0; length < encoded.size(); length++)
    {
        EXPECT_EQ(serializer.decode_spectra(const_cast<char*>(encoded.data()), length), nullptr)
            << "accepted prefix length " << length;
    }
}

TEST(BasicSerializer, RejectsOutOfRangeSpectrumIndex)
{
    Serializer serializer;
    std::string encoded = serializer.encode_spectra(make_block());

    // The final sparse pair is (index, value). Replace its index with size(),
    // which must be rejected instead of indexing beyond the Eigen array.
    const unsigned short invalid_index = 4;
    const size_t final_pair_size = sizeof(unsigned short) + sizeof(float);
    ASSERT_GE(encoded.size(), final_pair_size);
    std::memcpy(encoded.data() + encoded.size() - final_pair_size,
                &invalid_index, sizeof(invalid_index));

    EXPECT_EQ(serializer.decode_spectra(encoded.data(), encoded.size()), nullptr);
}

TEST(BasicSerializer, CountsRoundTripAndRejectsTruncation)
{
    Serializer serializer;
    const std::string encoded = serializer.encode_counts(make_block());
    auto decoded = serializer.decode_counts(const_cast<char*>(encoded.data()), encoded.size());

    ASSERT_NE(decoded, nullptr);
    ASSERT_TRUE(decoded->fitting_blocks.contains(data_struct::Fitting_Routines::ROI));
    EXPECT_FLOAT_EQ(decoded->fitting_blocks.at(data_struct::Fitting_Routines::ROI).fit_counts.at("Fe"), 42.0f);

    for (size_t length = 0; length < encoded.size(); length++)
    {
        EXPECT_EQ(serializer.decode_counts(const_cast<char*>(encoded.data()), length), nullptr)
            << "accepted prefix length " << length;
    }
}

TEST(BasicSerializer, RejectsUnknownFittingRoutine)
{
    Serializer serializer;
    std::string encoded = serializer.encode_counts(make_block());
    const size_t metadata_size = sizeof(unsigned int) + (4 * sizeof(size_t)) + sizeof(float) +
                                 std::string("scan.mda").size() + 1 +
                                 std::string("dataset/dir").size() + 1;
    const size_t routine_offset = metadata_size + sizeof(std::uint32_t);
    ASSERT_GE(encoded.size(), routine_offset + sizeof(std::uint32_t));

    const std::uint32_t invalid_routine = 0xffffffffU;
    std::memcpy(encoded.data() + routine_offset, &invalid_routine, sizeof(invalid_routine));
    EXPECT_EQ(serializer.decode_counts(encoded.data(), encoded.size()), nullptr);
}

TEST(BasicSerializer, CombinedDecoderRequiresSpectraSection)
{
    Serializer serializer;
    std::string counts_only = serializer.encode_counts(make_block());
    EXPECT_EQ(serializer.decode_counts_and_spectra(counts_only.data(), counts_only.size()), nullptr);
}

} // namespace
