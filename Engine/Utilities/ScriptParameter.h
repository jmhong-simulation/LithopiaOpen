#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>

class ScriptParameter
{
public:
	enum class Type { UNKNOWN, INT, FLOAT, BOOL, DOUBLE, STRING, FLOAT2, FLOAT3, FLOAT4, INT2, INT3, INT4};	// use UNKNOWN for argument reader	

	Type		type_;
	std::string	name_;
	std::string	value_;

	ScriptParameter()
	{}

	ScriptParameter(const Type _type, const std::string _name, const std::string _value)
		: type_(_type), name_(_name), value_(_value)
	{}
};