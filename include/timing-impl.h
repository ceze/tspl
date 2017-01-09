/*
 * Copyright (c) 2008-2011 Zhang Ming (M. Zhang), zmjerry@163.com
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 2 or any later version.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details. A copy of the GNU General Public License is available at:
 * http://www.fsf.org/licensing/licenses
 */


/*****************************************************************************
 *                               timing-impl.h
 *
 * Implementation for Timing class.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructor and destructor
 */
Timing::Timing() : running(false), startTime(0), total(0)
{
}

Timing::~Timing()
{
}


/**
 * timing start
 */
inline void Timing::start()
{
	running = true;
	total = 0;
	startTime = seconds();
}


/**
 * timing stop
 */
inline void Timing::stop()
{
	if( running )
	{
		 total += ( seconds() - startTime );
		 running = false;
	}
}


/**
 * Get the time between start() and stop().
 */
inline double Timing::read()
{
	if( running )
		stop();

	return total;
}
