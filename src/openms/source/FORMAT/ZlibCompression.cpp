// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ZlibCompression.h>

#include <QtCore/QByteArray>

#include <OpenMS/CONCEPT/LogStream.h>

#include <array>
#include <zlib.h>

using namespace std;

namespace OpenMS
{

  void ZlibCompression::compressString(std::string& str, std::string& compressed)
  {
    compressData(reinterpret_cast<Bytef*>(&str[0]), str.size(), compressed);
  }

  void ZlibCompression::compressData(const void* raw_data, const size_t in_length, std::string& compressed)
  {
    compressed.clear();

    const unsigned long sourceLen = (unsigned long)in_length;
    unsigned long compressed_length =                         // compressBound((unsigned long)str.size());
      sourceLen + (sourceLen >> 12) + (sourceLen >> 14) + 11; // taken from zlib's compress.c, as we cannot use compressBound*

    int zlib_error;

    compressed.resize(compressed_length); // reserve enough space -- we may not need all of it
    zlib_error = compress(reinterpret_cast<Bytef*>(&compressed[0]), &compressed_length, (Bytef*)raw_data, sourceLen);

    switch (zlib_error)
    {
      case Z_MEM_ERROR:
      case Z_BUF_ERROR:
        throw Exception::OutOfMemory(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, compressed_length);
      case Z_OK: // ok
        break;
    }

    if (zlib_error != Z_OK)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Compression error?");
    }
    compressed.resize(compressed_length); // cut down to the actual data
  }

  void ZlibCompression::compressString(const QByteArray& raw_data, QByteArray& compressed_data)
  {
    compressed_data = qCompress(raw_data);
    compressed_data.remove(0, 4);
  }

  void ZlibCompression::uncompressString(const void* compressed_data, size_t nr_bytes, std::string& raw_data, size_t output_size)
  {
    raw_data.resize(output_size);
    uLongf uncompressedSize = output_size;
    int ret = uncompress((Bytef*)raw_data.data(), &uncompressedSize, (Bytef*)compressed_data, nr_bytes);

    if (ret == Z_OK)
    {
      if (uncompressedSize != raw_data.size())
      { 
        OPENMS_LOG_INFO << "zlib::uncompress: data was smaller than anticipated.\n";
        raw_data.resize(output_size);
      }
    }
    else {
      std::cerr << "Zlib::uncompress() failed with code: " << ret << " and expected output size: " << output_size << std::endl;
    }

  }

  void ZlibCompression::uncompressString(const void* compressed_data, size_t nr_bytes, std::string& uncompressed)
  {
    const size_t CHUNK_SIZE = 16384;
    uncompressed.clear();
    z_stream strm = {};

    // Setup input
    strm.next_in = (Bytef*)(compressed_data);
    strm.avail_in = nr_bytes;

    // Initialize zlib (use inflateInit2 for gzip or raw deflate)
    if (inflateInit(&strm) != Z_OK) { throw std::runtime_error("inflateInit failed"); }

    // Decompress loop
    std::array<char, CHUNK_SIZE> buffer;
    int ret;

    do
    {
      strm.avail_out = CHUNK_SIZE;
      strm.next_out = (Bytef*)buffer.data();

      ret = inflate(&strm, Z_NO_FLUSH);
      if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR)
      {
        inflateEnd(&strm);
        throw std::runtime_error("inflate failed");
      }

      size_t bytesDecompressed = CHUNK_SIZE - strm.avail_out;
      uncompressed.insert(uncompressed.end(), buffer.begin(), buffer.begin() + bytesDecompressed);

    } while (ret != Z_STREAM_END);

    inflateEnd(&strm);
  }

  void ZlibCompression::uncompressString(const QByteArray& compressed_data, QByteArray& raw_data)
  {
    std::string uncompressed;
    uncompressString(compressed_data.constData(), compressed_data.size(), uncompressed);
    raw_data = QByteArray::fromRawData(uncompressed.data(), static_cast<int>(uncompressed.size()));
  }

}

