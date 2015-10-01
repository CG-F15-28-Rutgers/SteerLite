//
// Copyright (c) 2015 Mahyar Khayatkhoei
// Copyright (c) 2009-2014 Shawn Singh, Glen Berseth, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"

using namespace Util;

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{
#ifdef ENABLE_GUI

	//================DELETE THIS PART AND THEN START CODING===================
	static bool flag = false;
	if (!flag)
	{
		std::cerr << "ERROR>>>>Member function drawCurve is not implemented!" << std::endl;
		flag = true;
	}
	//=========================================================================

	// Robustness: make sure there is at least two control point: start and end points

	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points

	return;
#endif
}

// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	/* simple bubblesort control points */
	/* TODO(ahh42): optimize? */
	int len = controlPoints.size();
	bool swapped = false;
	if(len < 2)
	{
		return;
	}
	do
	{
		swapped = false;
		for (int i = 1; i < len; i++)
		{
			if (controlPoints[i].time < controlPoints[i - 1].time)
			{
				CurvePoint temp = controlPoints[i];
				controlPoints[i] = controlPoints[i - 1];
				controlPoints[i - 1] = temp;
				swapped = true;
			}
		}
	} while (swapped);
#ifdef DEBUG
	printf("Sorted Control Points: \n");
	for (int i = 0; i < len; i++)
	{
		//printf("%f %f %f\n", controlPoints[i].position.x, controlPoints[i].position.y, controlPoints[i].position.z);
		printf("\t%f\n", controlPoints[i].time);
	}
	printf("\n");
#endif

	return;
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;

	// Define temporary parameters for calculation
	unsigned int nextPoint;
	float normalTime, intervalTime;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Roboustness
bool Curve::checkRobust()
{
	if(controlPoints.size() < 2)
		return false;

	return true;
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	int len = controlPoints.size();
	if (time < 0.0f)
		return false;
	if(len < 1)
		return false;
	if(time < controlPoints[0].time)
	{
		nextPoint = 0;
		return true;
	}

	/* if the time is at or past the last control point return false */
	if (time >= controlPoints[len-1].time)
		return false;

	for (int i = 1; i < len; i++)
	{
		if (time < controlPoints[i].time && time >= controlPoints[i-1].time)
		{
			nextPoint = i;
#ifdef DEBUG
			printf("Current time: %f.\nChosen point index: %d, time: %f.\n", time, nextPoint, controlPoints[nextPoint].time);
/*			for (int i = 0; i < len; i++)
			{
				printf("\t%f\n", controlPoints[i].time);
			}
*/
#endif
			return true;
		}
	}

	return false;
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;

	// Calculate time interval, and normal time required for later curve calculations
	intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
	/*
	 * normalTime is the actual time normalized to be between 0 and 1, 0 being the time of the last control
	 * point and 1 being the time of the next control point
	*/
	normalTime = (time - controlPoints[nextPoint - 1].time) / (intervalTime);
#ifdef DEBUG
	printf("Normal time: %f.\n", normalTime);
#endif
	if (normalTime < 0 || normalTime > 1)
		return newPosition;

	Point p0 = controlPoints[nextPoint - 1].position;
	Point p1 = controlPoints[nextPoint].position;
	Vector t0 = controlPoints[nextPoint - 1].tangent * intervalTime;
	Vector t1 = controlPoints[nextPoint].tangent * intervalTime;
	// Calculate position at t = time on Hermite curve
	/* based on the formula in the lecture slides (curves pg 64)... */
	newPosition = (2*pow(normalTime, 3) - 3*pow(normalTime, 2) + 1)*p0 + (pow(normalTime, 3) - 2*pow(normalTime, 2) + normalTime)*t0 + (-2*pow(normalTime, 3) + 3*pow(normalTime, 2))*p1 + (pow(normalTime, 3) - pow(normalTime, 2))*t1;

#ifdef DEBUG
	printf("Moving toward control point %f %f %f.\n", controlPoints[nextPoint].position.x, controlPoints[nextPoint].position.y, controlPoints[nextPoint].position.z);
	printf("New position %f %f %f on curve.\n", newPosition.x, newPosition.y, newPosition.z);
#endif
	// Return result
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;

	//================DELETE THIS PART AND THEN START CODING===================
	static bool flag = false;
	if (!flag)
	{
		std::cerr << "ERROR>>>>Member function useCatmullCurve is not implemented!" << std::endl;
		flag = true;
	}
	//=========================================================================


	// Calculate time interval, and normal time required for later curve calculations

	// Calculate position at t = time on Catmull-Rom curve

	// Return result
	return newPosition;
}
