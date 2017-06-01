#pragma once

#include <stdio.h>
#include <boost/chrono.hpp>

class Timer
{
public:
	boost::chrono::system_clock::time_point start_time_;
	boost::chrono::duration<double> accumulated_time_;

public:
	Timer(void)
	{
		accumulated_time_ = boost::chrono::duration<double>(boost::chrono::duration_values<double>::zero());
	}

	~Timer(void)
	{}

public:
	void start()
	{
		printf("Timer start\n");
		start_time_ = boost::chrono::system_clock::now();
	}

	void end()
	{
		boost::chrono::duration<double> elapsed_time = boost::chrono::system_clock::now() - start_time_;

		printf("Timer end %f\n", (float)elapsed_time.count());

		accumulated_time_ += elapsed_time;
	}

	void coutAccumulatedTime()
	{
		printf("Accumulated time %f\n", (float)accumulated_time_.count());
	}
};



