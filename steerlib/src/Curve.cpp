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


	// Robustness: make sure there is at least two control point: start and end points
	checkRobust();

	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	float time = 0.0f + window;
	Point lastDrawPoint = controlPoints[0].position;
	while (time <= controlPoints[controlPoints.size()-1].time)
	{
		Point newDrawPoint;
		unsigned int nextPoint; 
		findTimeInterval(nextPoint, time);
		if (type == hermiteCurve)
		{
			newDrawPoint = useHermiteCurve(nextPoint, time);
		}
		else if (type == catmullCurve)
		{
			newDrawPoint = useCatmullCurve(nextPoint, time);
		}
		DrawLib::drawLine(lastDrawPoint, newDrawPoint, curveColor, curveThickness);
		time += window;
		//printf("time %f\n", time);
		lastDrawPoint = newDrawPoint;
	}

	return;
#endif
}

// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	/* simple bubblesort control points */
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
	/* remove control points that have the same time */
	for (int i = 0; i < controlPoints.size(); i++)
	{
		for (int j = i+1; j < controlPoints.size(); j++)
		{
			if (controlPoints[i].time == controlPoints[j].time)
				controlPoints.erase(controlPoints.begin() + j);
			else
				break;
		}
	}
#ifdef DEBUG
	printf("Sorted Control Points: \n");
	for (int i = 0; i < len; i++)
	{
		//printf("%f %f %f\n", controlPoints[i].position.x, controlPoints[i].position.y, controlPoints[i].position.z);
		//printf("\t%f\n", controlPoints[i].time);
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
	if (time < 0.0f)
		return false;
	if(controlPoints.size() < 1)
	{
		return false;
	}
	if(time < controlPoints[0].time)
	{
		return false;
	}

	/* if the time is at or past the last control point return false */
	if (time >= controlPoints[controlPoints.size()-1].time)
	{
		return false;
	}

	for (int i = 1; i <= controlPoints.size()-1; i++)
	{
		if (time < controlPoints[i].time && time >= controlPoints[i-1].time)
		{
			nextPoint = i;
#ifdef DEBUG
			printf("Current time: %f.\nChosen point index: %d, time: %f.\n", time, nextPoint, controlPoints[nextPoint].time);
			for (int i = 0; i < controlPoints.size(); i++)
			{
				//printf("\t%f %f %f %f\n", controlPoints[i].time, controlPoints[i].position.x, controlPoints[i].position.y, controlPoints[i].position.z);
			}

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

	Point p1 = controlPoints[nextPoint - 1].position;
	Point p2 = controlPoints[nextPoint].position;
	/* tangents (need to be scaled by time interval between last control point and next control point) */
	Vector s0 = controlPoints[nextPoint - 1].tangent * intervalTime;
	Vector s1 = controlPoints[nextPoint].tangent * intervalTime;
	/* based on the formula in the lecture slides (curves pg 64)... */
	newPosition = (2*pow(normalTime, 3) - 3*pow(normalTime, 2) + 1)*p1 + (pow(normalTime, 3) - 2*pow(normalTime, 2) + normalTime)*s0 + 
		(-2*pow(normalTime, 3) + 3*pow(normalTime, 2))*p2 + (pow(normalTime, 3) - pow(normalTime, 2))*s1;

#ifdef DEBUG
	printf("Moving toward control point %f %f %f.\n", controlPoints[nextPoint].position.x, controlPoints[nextPoint].position.y, controlPoints[nextPoint].position.z);
	//printf("New position %f %f %f on curve.\n", newPosition.x, newPosition.y, newPosition.z);
#endif
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	/* Catmull-Rom is just a hermite curve with a different tangent calculation (I think?) */
	Point newPosition;
	float normalTime, intervalTime;
	Point p1, p2;
	Vector s1, s2;

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
	/*tangent calculation on the first and last control points are different */
	if (nextPoint == 1)
	{
		p1 = controlPoints[nextPoint - 1].position;
		p2 = controlPoints[nextPoint].position;
		Point p3 = controlPoints[nextPoint + 1].position;
		s1 = (2*((p2 - p1) / intervalTime) - ((p3 - p1) / (2*intervalTime))) * intervalTime;
		s2 = ((p3 - p1) / (2*intervalTime)) * intervalTime;
	}
	else if (nextPoint == controlPoints.size()-1)
	{
		Point p0 = controlPoints[nextPoint - 2].position;
		p1 = controlPoints[nextPoint - 1].position;
		p2 = controlPoints[nextPoint].position;
		//Point p3 = controlPoints[nextPoint + 1].position;
		s1 = ((p2 - p0) / (2*intervalTime)) * intervalTime;
		//s2 = (2*((p2 - p1) / intervalTime) - ((p3 - p1) / (2*intervalTime))) * intervalTime;
		s2 = (2*((p2 - p0) / intervalTime) - ((p1 - p0) / (2*intervalTime))) * intervalTime;
	}
	else
	{
		Point p0 = controlPoints[nextPoint - 2].position;
		p1 = controlPoints[nextPoint - 1].position;
		p2 = controlPoints[nextPoint].position;
		Point p3 = controlPoints[nextPoint + 1].position;
		s1 = ((p2 - p0) / (2*intervalTime)) * intervalTime;
		s2 = ((p3 - p1) / (2*intervalTime)) * intervalTime;
	}
	newPosition = (2*pow(normalTime, 3) - 3*pow(normalTime, 2) + 1)*p1 + (pow(normalTime, 3) - 2*pow(normalTime, 2) + normalTime)*s1 + 
		(-2*pow(normalTime, 3) + 3*pow(normalTime, 2))*p2 + (pow(normalTime, 3) - pow(normalTime, 2))*s2;


	return newPosition;
}
