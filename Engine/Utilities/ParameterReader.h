#pragma once

#include <exprtk/exprtk.hpp>
#include "ArgParameterList.h"
#include "DataStructure/Vector2D.h"
#include "DataStructure/Vector3D.h"
#include "DataStructure/Vector4D.h"
#include "ScriptParameterList.h"

#define DEFINE_PARAMETER(type, name) const type name = par_reader.getValue(#name, type());
#define DEFINE_PARAMETER_WITH_DEFAULT(type, name, default_name) const type temp = par_reader.getValue(#name, type()); const type name = temp.compare(#default_name) == 0 ? default_name : temp;

class ParameterReader
{
public:
	ScriptParameterList argument_reader_;
	ScriptParameterList script_reader_;

	ParameterReader(const int argc, char *argv[], const bool _use_cout = false)
	{
		argument_reader_.initializeFromArguments(argc, argv);

		const std::string script_filename = argument_reader_.isValid(std::string("script")) ? argument_reader_.getStringValue(std::string("script")) : std::string("./image_arrangement.parameters");

		script_reader_.initialize(script_filename.c_str(), _use_cout);
	}

public:
	std::string getStringValue(const std::string par_name)
	{
		if (argument_reader_.isValid(std::string(par_name)))
		{
			return argument_reader_.getStringValue(std::string(par_name));
		}
		else
		{
			return script_reader_.getStringValue(par_name);
		}
	}

	bool getValue(const std::string par_name, const bool type_specifier)
	{
		std::string value = getStringValue(par_name);

		if (value.compare("false") == 0 || value.compare("0") == 0) return false;
		else return true;
	}

	int getValue(const std::string par_name, const int type_specifier)
	{
		return (int)getValueFromString(getStringValue(par_name));
	}

	std::string getValue(const std::string par_name, const std::string type_specifier)
	{
		return getStringValue(par_name);
	}

	float getValue(const std::string par_name, const float type_specifier)
	{
		return getValueFromString(getStringValue(par_name));
	}

	Vector2D<float> getValue(const std::string par_name, const Vector2D<float> type_specifier)
	{
		std::vector<std::string> value_splits;

		std::string value = getStringValue(par_name);

		boost::split(value_splits, value, boost::is_any_of(" "), boost::token_compress_on);

		return Vector2D<float>(getValueFromString(value_splits[0]), getValueFromString(value_splits[1]));
	}

	Vector4D<float> getValue(const std::string par_name, const Vector4D<float> type_specifier)
	{
		std::vector<std::string> value_splits;

		std::string value = getStringValue(par_name);

		boost::split(value_splits, value, boost::is_any_of(" "), boost::token_compress_on);

		return Vector4D<float>(getValueFromString(value_splits[0]), getValueFromString(value_splits[1]), getValueFromString(value_splits[2]), getValueFromString(value_splits[3]));
	}

	float getValueFromString(const std::string string_expression) const
	{
		exprtk::expression<float> expression;
		exprtk::parser<float>     parser;

		parser.compile(string_expression, expression);

		return expression.value();
	}
};