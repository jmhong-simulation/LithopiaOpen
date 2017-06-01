// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <iostream>
#include <fstream>
#include <list>
#include "ScriptParameter.h"
#include "DataStructure/Vector2D.h"
#include "DataStructure/Vector3D.h"
#include "DataStructure/Vector4D.h"

#define DEFINE_PARAMETER(type, name) const type name = par_reader.getValue(#name, type());
#define DEFINE_PARAMETER_WITH_DEFAULT(type, name, default_name) const type temp = par_reader.getValue(#name, type()); const type name = temp.compare(#default_name) == 0 ? default_name : temp;

class ScriptParameterList
{
public:
	bool use_cout_;

	std::list<ScriptParameter> par_list_;

	ScriptParameterList();
	ScriptParameterList(const char* filename, const bool _use_cout = false);

	void		initialize(const char* filename, const bool _use_cout = false);
	void		initializeFromArgumentsAndFile(const int argc, char* argv[], const std::string default_script_filename);
	void		initializeFromArguments(const int argc, char* argv[]);
	void		initializeFromMessage(const std::string message);
	void		initializeFromFile(const char *filename);
	void		initializeBlockFromFile(const char *filename, std::vector<std::string>& block);
	void		overwriteFromArtuments(const int argc, char* argv[]);
	void		combineWith(ScriptParameterList& input_list);
	void		parseBlock(const std::vector<std::string>& block);
	void		parseLine(const std::string& line);
	void		parseTokens(std::list<std::string>& tokens);
	void		printAll();

	bool		isValid(const std::string par_name);	

	std::string getStringValue(const std::string par_name);
	void		setStringValue(const std::string par_name, const std::string par_value);

	int				getValue(const std::string par_name, const int type_specifier);
	bool			getValue(const std::string par_name, const bool type_specifier);
	float			getValue(const std::string par_name, const float type_specifier);
	std::string		getValue(const std::string par_name, const std::string type_specifier);
	Vector2D<int> getValue(const std::string par_name, const Vector2D<int> type_specifier);
	Vector3D<int> getValue(const std::string par_name, const Vector3D<int> type_specifier);
	Vector4D<int> getValue(const std::string par_name, const Vector4D<int> type_specifier);
	Vector2D<float> getValue(const std::string par_name, const Vector2D<float> type_specifier);
	Vector3D<float> getValue(const std::string par_name, const Vector3D<float> type_specifier);
	Vector4D<float> getValue(const std::string par_name, const Vector4D<float> type_specifier);

	float getValueFromString(const std::string string_expression) const;
};