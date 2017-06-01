// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "ScriptParameterList.h"

class ScriptBlock
{
public:
	std::string name_;

	ScriptParameterList par_list_;

	ScriptBlock *parent_;
	std::list<ScriptBlock*> children_;		// list of sub blocks

	ScriptBlock();
	~ScriptBlock();

	void parseTokens(std::list<std::string>& tokens_);

	ScriptBlock& getBlock(const std::string block_name)
	{
		for (std::list<ScriptBlock*>::iterator itr = children_.begin(); itr != children_.end(); ++itr)
		{
			if (block_name.compare((*itr)->name_) == 0) return (*(*itr));
		}

		std::cout << "block " << block_name << " was not found" << std::endl;

		exit(1);

		return *this;		
	}

	void print()
	{
		std::cout << "Block "<< name_ << std::endl;
		par_list_.printAll();

		for (std::list<ScriptBlock*>::iterator itr = children_.begin(); itr != children_.end(); ++itr)
			(*itr)->print();
	}

	template<class TT>
	TT getValue(const std::string& par_name, const TT& type_specifier)
	{
		return par_list_.getValue(par_name, type_specifier);
	}
};