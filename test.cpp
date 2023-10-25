#include "DetectorMovementArbitration.h"
#include <iostream>

int main()
{
	MovementArbitration ma;
	//Point detector(-5.5753, 10.0581),
	//	approach_point(-4, 2),
	//	start_point(-7.2296, 6.8545),
	//	end_point(-1.9818, 9.7633);
	//Point result_point = ma.RadiusDirectionApproach(detector, approach_point, start_point, end_point);
	//std::cout << result_point.x << " " << result_point.y << std::endl;
	Point vec_detect(-3.3452, 6.0349);

	std::vector<Point> convex_a =
	{
		{-7.3769,6.5511},
		{-4.9142,3.3990},
		{-0.1862,7.093},
		{-2.6488,10.245}
	};

	std::vector<Point> convex_b =
	{
		{-4,2},
		{-5.0934,-0.7755},
		{-3.7662,-3.4571},
		{-1.4367,-5.1906},
		{5.8224,-4.703},
		{5.9037,-0.4776},
		{4,4},
		{-0.9221,3.1791}
	};

	//ma.MinkowskiDifferenceSpace(convex_a, convex_b);
	//double d = ma.QueryTheMinkowskiSpace(convex_a, convex_b, 20 * M_PI / 180);
	//std::cout << d << std::endl;
	 
	//auto d = ma.MainRotateApproach(convex_a, convex_b);
	//std::cout << std::get<0>(d) << " " << std::get<1>(d) << std::endl;

	auto d = ma.Spin2Approach(convex_a, convex_b);
	std::cout << std::get<0>(d) << " " << std::get<1>(d) << std::endl;
	return 0;
}