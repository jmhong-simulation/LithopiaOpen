#include "ArgParameterList.h"
#include <boost/algorithm/string.hpp>

ArgParameterList::ArgParameterList(const int argc, char *argv[])
{
	readFromArguments(argc, argv);
}

void ArgParameterList::readFromArguments(const int argc, char *argv[])
{
	const std::string exeutable = argv[0];		// name of execution (not used)

	arg_list_.reserve((argc - 1) / 2);

	if ((argc - 1) % 2 != 0)	// argc should be a odd number. arguments are : name of execution -flag1 value1 -flag2 value2 ...
	{
		std::cout << "Wrong arguments" << std::endl;
		std::cin.get();
		exit(1);
	}

	// read flag (-name) and value pairs
	int i = 1;
	while (i < argc)
	{
		std::string flag = argv[i++];
		const std::string value = argv[i++];

		if (!boost::istarts_with(flag, "-"))
		{
			std::cout << "Flag " << flag << " does not start with -" << std::endl;
			exit(1);
		}

		boost::replace_first(flag, "-", "");	// remove "-" from flag to obtain parameter name. flag is "-name".

		arg_list_.push_back(ScriptParameter(ScriptParameter::Type::UNKNOWN, flag, value));
	}
}

std::string ArgParameterList::getStringValue(const std::string name) const
{
	for (int i = 0; i < arg_list_.size(); ++i)
	{
		if (arg_list_[i].name_.compare(name) == 0) return arg_list_[i].value_;
	}

	std::cout << "Argument " << name << " was not found. Terminating." << std::endl;
	std::cin.get();
	exit(1);
}

bool ArgParameterList::isValid(const std::string flag) const
{
	for (int i = 0; i < arg_list_.size(); ++i)
	{
		if (arg_list_[i].name_.compare(flag) == 0) return true;
	}

	return false;
}