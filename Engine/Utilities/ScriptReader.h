// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <list>
#include <fstream>
#include "ScriptBlock.h"

class ScriptReader
{
public:
	std::list<std::string> keywords_;
	std::list<std::string> tokens_;

	ScriptBlock *head_block_;

	ScriptReader()
		: head_block_(nullptr)
	{
		keywords_.push_back("int");
		keywords_.push_back("int2");
		keywords_.push_back("int3");
		keywords_.push_back("float");
		keywords_.push_back("float2");
		keywords_.push_back("float3");
		keywords_.push_back("float4");
		keywords_.push_back("bool");
		keywords_.push_back("double");
		keywords_.push_back("string");		
	}

	~ScriptReader()
	{
		SAFE_DELETE(head_block_);
	}

	void readFile(const char* filename);
	void findAllTokens(std::ifstream& file);
	void parseTokens();
	void parseTokens(std::list<std::string>& tokens_, ScriptBlock* block);
	bool isKeyword(const std::string& input);
	void printAll();
	void printBlock(ScriptBlock* block);

	ScriptBlock& getBlock(const std::string block_name)
	{
		return head_block_->getBlock(block_name);
	}

	//TODO: delete tokens and keywords after reading
};
