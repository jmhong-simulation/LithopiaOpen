// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "IMAGE_2D.h"
#include "../DataStructure/GridUniform2D.h"
#include <IL/il.h>
#include <IL/ilu.h>
#include <algorithm>

const bool IMAGE_2D::WriteBMP24(const char *filename)
{
	int i, j, ipos;
	int bytesPerLine;
	unsigned char *buffer;

	FILE *file;
	struct BMPHeader header;
	
	int width = res_x_;
	int height = res_y_;

	// The length of each line must be a multiple of 4 bytes
	bytesPerLine = (3 * (width + 1) / 4) * 4;

	strcpy(header.file_type_, "BM");
	header.offBits_ = 54;
	header.file_size_ = header.offBits_ + bytesPerLine * height;
	header.reserved_ = 0;
	header.header_size_ = 40;
	header.width_ = width;
	header.height_ = height;
	header.planes_ = 1;
	header.bit_count_ = 24;
	header.compression_ = 0;
	header.size_image_ = bytesPerLine * height;
	header.pixel_per_meter_x_ = 0;
	header.pixel_per_meter_y_ = 0;
	header.num_color_palette_ = 0;
	header.num_color_important_ = 0;

	file = fopen(filename, "wb");
	if (file == NULL) return(0);

	fwrite(&header.file_type_, 2, 1, file);
	fwrite(&header.file_size_, 4, 1, file);
	fwrite(&header.reserved_, 4, 1, file);
	fwrite(&header.offBits_, 4, 1, file);
	fwrite(&header.header_size_, 4, 1, file);
	fwrite(&header.width_, 4, 1, file);
	fwrite(&header.height_, 4, 1, file);
	fwrite(&header.planes_, 2, 1, file);
	fwrite(&header.bit_count_, 2, 1, file);
	fwrite(&header.compression_, 4, 1, file);
	fwrite(&header.size_image_, 4, 1, file);
	fwrite(&header.pixel_per_meter_x_, 4, 1, file);
	fwrite(&header.pixel_per_meter_y_, 4, 1, file);
	fwrite(&header.num_color_palette_, 4, 1, file);
	fwrite(&header.num_color_important_, 4, 1, file);

	buffer = (unsigned char*)malloc(bytesPerLine);
	if (buffer == NULL)
	{
		fprintf(stderr, "Can't allocate memory for BMP file.\n");
		return false;
	}

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			ipos = 3 * (width * i + j);
			buffer[3 * j] = data_(j, i).b_;
			buffer[3 * j + 1] = data_(j, i).g_;
			buffer[3 * j + 2] = data_(j, i).r_;
		}
		fwrite(buffer, bytesPerLine, 1, file);
	}

	free(buffer);
	fclose(file);

	return true;
}

const bool IMAGE_2D::ReadBMP24(const char * imagepath)
{
	printf("Reading image %s\n", imagepath);

	// Data read from the header of the BMP file
	unsigned char header[54];
	unsigned int dataPos;
	unsigned int imageSize;
	unsigned int width, height;

	// Open the file
	FILE * file = fopen(imagepath, "rb");
	if (!file) { printf("%s could not be opened. Are you in the right directory ? Don't forget to read the FAQ !\n", imagepath); getchar(); return 0; }

	// Read the header, i.e. the 54 first bytes

	// If less than 54 bytes are read, problem
	if (fread(header, 1, 54, file) != 54){
		printf("Not a correct BMP file\n");
		return false;
	}
	// A BMP files always begins with "BM"
	if (header[0] != 'B' || header[1] != 'M'){
		printf("Not a correct BMP file\n");
		return false;
	}

	// Make sure this is a 24bpp file
	if (*(int*)&(header[0x1E]) != 0)         { printf("Not a correct BMP file\n");    return false; }
	if (*(int*)&(header[0x1C]) != 24)         { printf("Not a correct BMP file\n");    return false; }

	// Read the information about the image
	dataPos = *(int*)&(header[0x0A]);
	imageSize = *(int*)&(header[0x22]);
	width = *(int*)&(header[0x12]);
	height = *(int*)&(header[0x16]);

	// Some BMP files are mis-formatted, guess missing information
//	if (imageSize == 0)    imageSize = width*height*3; // 3 : one byte for each Red, Green and Blue component
	if (dataPos == 0)      dataPos = 54; // The BMP header is done that way

	res_x_  = width;
	res_y_ = height;

	int scanline_byte = res_x_ * 3;
	int padding = 0;

	while((scanline_byte+padding) % 4 != 0) padding++;

	int psb = scanline_byte + padding; // padded scanline byte
	imageSize = psb*height;

	// Read the actual data from the file into the buffer
	unsigned char* image_buf = new unsigned char[imageSize];
	fread(image_buf, 1, imageSize, file);

	// Create a buffer
	data_.initialize(0, 0, width, height, false);

	long buf_pos = 0;
	long new_pos = 0;

	for(int y=0; y<height; y++ )
	for(int x=0; x<3*width; x+=3 )
	{
		new_pos = y*width + x/3;     
		buf_pos = y*psb + x;

		data_.values_[new_pos].r_ = image_buf[buf_pos+2];       
		data_.values_[new_pos].g_ = image_buf[buf_pos+1]; 
		data_.values_[new_pos].b_ = image_buf[buf_pos+0];     
	}

	// Everything is in memory now, the file wan be closed
	fclose(file);

	delete image_buf;

	return true;
}

const int IMAGE_2D::GetResizeWidth(const int new_height) const
{
	const T rho = (T)new_height / (T)res_y_;

	return (int)floor((T)res_x_ * rho);
}

const int IMAGE_2D::GetResizeHeight(const int new_width) const
{
	const T rho = (T)new_width / (T)res_x_;

	return (int)floor((T)res_y_ * rho);
}

void IMAGE_2D::ScaleTo(const int new_i_res, const int new_j_res)
{
	GridUniform2D grid_original;
	grid_original.Initialize(data_.i_start_, data_.j_start_, data_.i_res_, data_.j_res_, 0, 0, 1, 1);

	GridUniform2D grid_scaled;
	grid_scaled.Initialize(0, 0, new_i_res, new_j_res, 0, 0, 1, 1);

	Array2D<BMP_RGB> array_temp(data_);

	data_.initialize(0, 0, new_i_res, new_j_res, 0);

	for (int j = data_.j_start_; j <= data_.j_end_; ++j)
	for (int i = data_.i_start_; i <= data_.i_end_; ++i)
	{
		data_(i, j) = grid_original.GetClampedLinearInterpolationCell(array_temp, grid_scaled.GetCellCenter(i, j));
	}

	res_x_ = new_i_res;
	res_y_ = new_j_res;
}

void IMAGE_2D::ScaleBy(const T scale_x, const T scale_y)
{
	const int new_i_res = (int)floor((T)res_x_*scale_x);
	const int new_j_res = (int)floor((T)res_y_*scale_y);

	ScaleTo(new_i_res, new_j_res);
}

const T IMAGE_2D::GetAspectRatio() const
{
	return (T)width_ / (T)height_;
}

void IMAGE_2D::ScaleToFitHeight(const int new_height)
{
	ScaleTo(GetResizeWidth(new_height), new_height);
}

void IMAGE_2D::ScaleToFitWidth(const int new_width)
{
	ScaleTo(GetResizeHeight(new_width), new_width);
}

void IMAGE_2D::ScaleIso(const T scale)
{ 
	ScaleBy(scale, scale); 
}

void IMAGE_2D::Initialize(const int res_x_input, const int res_y_input)
{
	res_x_ = res_x_input;
	res_y_ = res_y_input;

	data_.initialize(0, 0, res_x_, res_y_, false);
}

void IMAGE_2D::Initialize(const IMAGE_2D& image)
{
	Initialize(image.width_, image.height_);

	data_.copyFrom(image.data_);
}

// initialize a grey scale image from a height map
void IMAGE_2D::Initialize(const Array2D<T>& height_map)
{
	Initialize(height_map.i_res_, height_map.j_res_);

	for(int j = 0; j < res_y_; ++j)
		for (int i = 0; i < res_x_; ++i)
			data_(i, j).SetGray(height_map(i, j));
}

void IMAGE_2D::Rotate90()
{
	Array2D<BMP_RGB> data_temp(0, 0, res_y_, res_x_, false);

	for (int j = 0; j < res_y_; ++j)
	for (int i = 0; i < res_x_; ++i)
	{
		data_temp(j, i).CopyFrom(data_(i, j));

//		std::cout << data_temp(j, i).GetReversedGrey() << std::endl;
	}

	SWAP(data_.values_, data_temp.values_, BMP_RGB*);

	int t = res_x_;
	res_x_ = res_y_;
	res_y_ = t;

	data_.swapIJ();
}

void IMAGE_2D::ReflectLeftRight()
{
	Array2D<BMP_RGB> data_temp(0, 0, res_x_, res_y_, false);

	for (int j = 0; j < res_y_; ++j)
	for (int i = 0; i < res_x_; ++i)
	{
		data_temp(res_x_ -1 - i, j).CopyFrom(data_(i, j));
	}

	SWAP(data_.values_, data_temp.values_, BMP_RGB*);
}

void IMAGE_2D::ReflectUpDown()
{
	Array2D<BMP_RGB> data_temp(0, 0, res_x_, res_y_, false);

	for (int j = 0; j < res_y_; ++j)
	for (int i = 0; i < res_x_; ++i)
	{
		data_temp(i, res_y_ - 1 - j).CopyFrom(data_(i, j));
	}

	SWAP(data_.values_, data_temp.values_, BMP_RGB*);
}

void IMAGE_2D::ReverseColors()
{
	for (int j = 0; j < res_y_; ++j)
	for (int i = 0; i < res_x_; ++i)
	{
		data_(i,j).ReverseColor();
	}
}

void IMAGE_2D::ExtendOneColumn(const int column, const int width)
{
	data_.extendOneColumn(column, width);

	res_x_ = data_.i_res_;
	res_y_ = data_.j_res_;
}

void IMAGE_2D::CutLeftAndRight(const int cut_width)
{
	assert(cut_width <= width_);

	const int cut_width_left = cut_width/2;

	Array2D<BMP_RGB> data_temp(0, 0, width_ - cut_width, height_, false);

	for (int j = 0; j < data_temp.j_res_; ++j)
	for (int i = 0; i < data_temp.i_res_; ++i)
	{
		data_temp(i, j).CopyFrom(data_(i + cut_width_left, j));
	}

	SWAP(data_.values_, data_temp.values_, BMP_RGB*);

	data_.i_res_ = data_temp.i_res_;
	data_.j_res_ = data_temp.j_res_;

	width_ = data_.i_res_;
	height_ = data_.j_res_;
}

void IMAGE_2D::CutTopAndBottom(const int cut_width)
{
	assert(cut_width <= width_);

	const int cut_width_bottom = cut_width / 2;

	Array2D<BMP_RGB> data_temp(0, 0, width_, height_ - cut_width, false);

	for (int j = 0; j < data_temp.j_res_; ++j)
	for (int i = 0; i < data_temp.i_res_; ++i)
	{
		data_temp(i, j).CopyFrom(data_(i, j + cut_width_bottom));
	}

	SWAP(data_.values_, data_temp.values_, BMP_RGB*);

	data_.i_res_ = data_temp.i_res_;
	data_.j_res_ = data_temp.j_res_;

	width_ = data_.i_res_;
	height_ = data_.j_res_;
}

const BOX_2D<int> IMAGE_2D::GetColorBox(const BMP_RGB& color_mask) const
{
	BOX_2D<int> box(INT_MAX, INT_MAX, INT_MIN, INT_MIN);

	for (int j = data_.j_start_; j <= data_.j_end_; ++j)
	for (int i = data_.i_start_; i <= data_.i_end_; ++i)
	{
		if (data_(i, j).IsEqual(color_mask)) box.Include(Vector2D<int>(i, j));
	}

	return box;
}

void IMAGE_2D::SetAllColor(const BMP_RGB& input_color)
{
	for (int j = 0; j < res_y_; ++j)
	for (int i = 0; i < res_x_; ++i)
	{
		data_(i, j) = input_color;
	}
}

void IMAGE_2D::SetMultipleAllColor(T val)
{
	for (int j = 0; j < res_y_; ++j)
	for (int i = 0; i < res_x_; ++i)
	{
		data_(i, j).r_ = (unsigned char)(val * (int)data_(i, j).r_);
		data_(i, j).g_ = (unsigned char)(val * (int)data_(i, j).g_);
		data_(i, j).b_ = (unsigned char)(val * (int)data_(i, j).b_);
	}
}

void IMAGE_2D::ReadFileAndFitHeight(const std::string image_file_path, const int target_height)
{
	ilInit();

	using namespace std;

	string image_file_name;	//TODO: remove

	// Find file name from path
	if ((int)image_file_path.find("\\") > -1 || (int)image_file_path.find("/") > -1)
	{
		int ix = max((int)image_file_path.find_last_of("\\"), (int)image_file_path.find_last_of("/"));	//TODO: need to include <algorithm> only for this?
		image_file_name = image_file_path.substr(ix + 1);
	}
	else
	{
		image_file_name = image_file_path;
	}

	if (ilLoadImage(image_file_path.c_str()) == false)
	{
		cout << "Cannot read image file " << image_file_path.c_str() << endl;
		exit(1);
	}

	const int input_width = ilGetInteger(IL_IMAGE_WIDTH);
	const int input_height = ilGetInteger(IL_IMAGE_HEIGHT);
	const int input_size = ilGetInteger(IL_IMAGE_SIZE_OF_DATA);
	const int input_ch = ilGetInteger(IL_IMAGE_CHANNELS);

	if (input_ch < 3)
	{
		cout << "Image channel < 3" << endl;
		exit(1);
	}

	const int ilu_target_height = target_height;
	const float rho = (float)input_height / (float)ilu_target_height;
	const int ilu_target_width = input_width / rho;

	iluImageParameter(ILU_FILTER, ILU_BILINEAR);
	iluScale(ilu_target_width, ilu_target_height, ilGetInteger(IL_IMAGE_DEPTH));

	Initialize(ilu_target_width, ilu_target_height);

	const ILubyte *bytes = ilGetData();

	const int image_width_temp = ilGetInteger(IL_IMAGE_WIDTH);
	const int image_height_temp = ilGetInteger(IL_IMAGE_HEIGHT);

	if (IL_BMP == ilTypeFromExt(image_file_name.c_str()))
	{
		// for bmp
		for (int j = 0; j < res_y_; j++)
		{
			for (int i = 0; i < res_x_; i++)
			{
				int ix = (j * res_x_ + i)*input_ch;

				data_(i, j).r_ = bytes[ix + 2];
				data_(i, j).g_ = bytes[ix + 1];
				data_(i, j).b_ = bytes[ix + 0];
			}
		}
	}
	else
	{
		// for others
		for (int j = 0; j < res_y_; j++)
		{
			for (int i = 0; i < res_x_; i++)
			{
				int ix = ((res_y_ - 1 - j)*res_x_ + i)*input_ch;

				data_(i, j).r_ = bytes[ix + 0];
				data_(i, j).g_ = bytes[ix + 1];
				data_(i, j).b_ = bytes[ix + 2];
			}
		}
	}

	//TODO: empty il memory after reading
}
