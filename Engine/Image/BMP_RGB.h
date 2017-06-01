// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "../Engine/GENERIC_DEFINITIONS.h"

struct BMPHeader
{
	// http://en.wikipedia.org/wiki/BMP_file_format

	char file_type_[2];        // "BM"
	int file_size_;            // Size of file in bytes
	int reserved_;             // set to 0
	int offBits_;              // Byte offset to actual bitmap data (= 54)
	int header_size_;          // Size of BITMAPINFOHEADER, in bytes (= 40)
	int width_;                // Width of image, in pixels
	int height_;               // Height of images, in pixels
	short planes_;             // Number of planes in target device (set to 1)
	short bit_count_;          // Bits per pixel (24 in this case)
	int compression_;          // Type of compression (0 if no compression)
	int size_image_;           // Image size, in bytes (0 if no compression)
	int pixel_per_meter_x_;    // Horizontal resolution of the image (pixel per meter)
	int pixel_per_meter_y_;    // Vertical resolution of the image (pixel per meter)
	int num_color_palette_;    // Number of colors in the color table (if 0, use maximum allowed by biBitCount)
	int num_color_important_;  // Number of important colors.  If 0, all colors are important
};

class BMP_RGB
{
public:
	unsigned char r_, g_, b_;

	BMP_RGB()
	{}

	BMP_RGB(const unsigned char r, const unsigned char g, const unsigned char b)
		: r_(r), g_(g), b_(b)
	{}

	void CopyFrom(const BMP_RGB& rgb)
	{
		r_ = rgb.r_;
		g_ = rgb.g_;
		b_ = rgb.b_;
	}

	const T GetGray() const
	{
		// luminosity http://www.johndcook.com/blog/2009/08/24/algorithms-convert-color-grayscale/

		static const T normalizer = 1.0f / 255.0f;

		return (0.21f*r_ + 0.72f*g_ + 0.07f*b_) * normalizer;
	}

	const void SetGray(const T& grey_color)
	{
		const unsigned char value = (int)(grey_color*255.0f);

		r_ = g_ = b_ = value;
	}

	const T GetReversedGray() const
	{
		return 1.0f - GetGray();
	}

	const void ReverseColor()
	{
		r_ = (unsigned char)((int)255 - (int)r_);
		g_ = (unsigned char)((int)255 - (int)g_);
		b_ = (unsigned char)((int)255 - (int)b_);
	}

	const bool IsBlue() const
	{
		if ((int)r_ == (int)0 && (int)g_ == (int)0 && (int)b_ == (int)255) return true;
		else return false;
	}

	const bool IsGreen() const
	{
		if ((int)r_ == (int)0 && (int)g_ == (int)255 && (int)b_ == (int)0) return true;
		else return false;
	}

	const bool IsBlack() const
	{
		if ((int)r_ == (int)0 && (int)g_ == (int)0 && (int)b_ == (int)0) return true;
		else return false;	
	}

	const bool IsWhite() const
	{
		if ((int)r_ == (int)255 && (int)g_ == (int)255 && (int)b_ == (int)255) return true;
		else return false;
	}

	const bool IsRed() const
	{
		if ((int)r_ == (int)255 && (int)g_ == (int)0 && (int)b_ == (int)0) return true;
		else return false;
	}

	const bool IsEqual(const BMP_RGB& rgb) const
	{
		if (r_ == rgb.r_ && g_ == rgb.g_ && b_ == rgb.b_) return true;
		else return false;
	}

	void SetWhite()
	{
		r_ = (unsigned char)255;
		g_ = (unsigned char)255;
		b_ = (unsigned char)255;
	}

	void SetRed()
	{
		r_ = (unsigned char)255;
		g_ = (unsigned char)0;
		b_ = (unsigned char)0;
	}

	void SetBlue()
	{
		r_ = (unsigned char)0;
		g_ = (unsigned char)0;
		b_ = (unsigned char)255;
	}

	void SetGreen()
	{
		r_ = (unsigned char)0;
		g_ = (unsigned char)255;
		b_ = (unsigned char)0;
	}
};

template<class TT> static BMP_RGB operator * (const TT& s, const BMP_RGB& bmp)
{
	return BMP_RGB(s*(TT)bmp.r_, s*(TT)bmp.g_, s*(TT)bmp.b_);
}

template<class TT> static BMP_RGB operator * (const BMP_RGB& bmp, const TT& s)
{
	return BMP_RGB(s*(TT)bmp.r_, s*(TT)bmp.g_, s*(TT)bmp.b_);
}

static BMP_RGB operator + (const BMP_RGB& bmp1, const BMP_RGB& bmp2)
{
	return BMP_RGB(bmp1.r_ + bmp2.r_, bmp1.g_ + bmp2.g_, bmp1.b_ + bmp2.b_);
}