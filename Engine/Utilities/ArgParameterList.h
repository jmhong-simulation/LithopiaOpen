#pragma once

#include <fstream>
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include "ScriptParameter.h"

class ArgParameterList
{
public:
	std::vector<ScriptParameter> arg_list_;

	ArgParameterList(const int argc, char *argv[]);
	
	bool		isValid(const std::string flag) const;
	void		readFromArguments(const int argc, char *argv[]);
	std::string getStringValue(const std::string flag) const;
};