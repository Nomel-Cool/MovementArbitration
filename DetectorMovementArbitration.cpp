#include "DetectorMovementArbitration.h"

Point MovementArbitration::RadiusDirectionApproach(const Point& vec_detector, const Point& approach_point, const Point& start_point, const Point& end_point, double safe_distance = 1)
{
	/* 计算边向量 */
	Point edge(end_point.x - start_point.x, end_point.y - start_point.y);

	/* 作点到边的向量 */
	Point vec_point_2_edge(approach_point.x - start_point.x, approach_point.y - start_point.y);

	/* 求上述两条向量的长度 */
	double length_edge = sqrt(edge.x * edge.x + edge.y * edge.y);
	double length_point2edge = sqrt(vec_point_2_edge.x * vec_point_2_edge.x + vec_point_2_edge.y * vec_point_2_edge.y);

	/* 作边向量到点边向量的叉乘 */
	double cross_result = Cross(edge, vec_point_2_edge);

	/* 求上述两向量的夹角sin值 */
	double sin_value = cross_result / (length_edge * length_point2edge);

	/* 求距离向量 */
	Point vec_distance(vec_point_2_edge.x * sin_value, vec_point_2_edge.y * sin_value);

	/* 求出可接近距离 */
	double distance = sqrt(vec_distance.x * vec_distance.x + vec_distance.y * vec_distance.y);
	double movable_distance = distance - safe_distance;

	/* 计算探头向量的长度 */
	double length_detector = sqrt(vec_detector.x * vec_detector.x + vec_detector.y * vec_detector.y);

	/* 求出距离向量与探头向量的cos值 */
	double cos_value = ProjectOn(vec_detector, vec_distance) / (distance * length_detector);

	/* 求出探头可移动距离 */
	double detector_movable_distance = abs(movable_distance * cos_value);

	/* 求出探头径向移动后的坐标 */
	Point vec_moved_detector(vec_detector.x - detector_movable_distance * vec_detector.x / length_detector, vec_detector.y - detector_movable_distance * vec_detector.y / length_detector);

	return vec_moved_detector;
}

std::tuple<double, double> MovementArbitration::Spin2Approach(const std::vector<Point>& detector, const std::vector<Point>& convex)
{
	Point origin(0, 0);
	Point min_vec = closest_point_to_origin(detector, origin);
	double min_len = sqrt(min_vec.x * min_vec.x + min_vec.y * min_vec.y);

	/* 预处理,作明可夫斯基差空间 */
	MinkowskiDifferenceSpace_Spin(detector, convex, 0);

	/* 步进查询，步进值为1度 */
	double left_max = FindMaxCCWAngle(detector, convex, M_PI * 0.01 / 180);
	double right_max = FindMaxCWAngle(detector, convex, M_PI * 0.01 / 180);

	return { left_max,right_max };
}

void MovementArbitration::MinkowskiDifferenceSpace_Rotate(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b)
{
	combine4cos.resize(convex_a.size());
	for (int g = 0; g < combine4cos.size(); g++)
		combine4cos[g].resize(convex_b.size());

	combine4sin.resize(convex_a.size());
	for (int g = 0; g < combine4sin.size(); g++)
		combine4sin[g].resize(convex_b.size());

	for (int i = 0; i < convex_a.size(); i++)
	{
		double length_a2 = convex_a[i].x * convex_a[i].x + convex_a[i].y * convex_a[i].y;
		list_len_a.push_back(length_a2);
	}
	for (int j = 0; j < convex_b.size(); j++)
	{
		double length_b2 = convex_b[j].x * convex_b[j].x + convex_b[j].y * convex_b[j].y;
		list_len_b.push_back(length_b2);
	}
	for (int i = 0; i < convex_a.size(); i++)
	{
		for (int j = 0; j < convex_b.size(); j++)
		{
			combine4cos[i][j] = -2 * (convex_a[i].x * convex_b[j].x + convex_a[i].y * convex_b[j].y);
			combine4sin[i][j] = 2 * (convex_a[i].y * convex_b[j].x - convex_a[i].x * convex_b[j].y);
		}
	}
}

void MovementArbitration::MinkowskiDifferenceSpace_Spin(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double theta)
{
	spin_len_a.resize(convex_a.size());

	spin_len_b.resize(convex_b.size());

	common_value.resize(convex_a.size());
	for (int g = 0; g < common_value.size(); g++)
		common_value[g].resize(convex_b.size());

	for (int i = 0; i < convex_a.size(); i++)
	{
		double length_a2 = convex_a[i].x * convex_a[i].x + convex_a[i].y * convex_a[i].y;
		spin_len_a[i] = length_a2;
	}
	for (int j = 0; j < convex_b.size(); j++)
	{
		double length_b2 = convex_b[j].x * convex_b[j].x + convex_b[j].y * convex_b[j].y;
		spin_len_b[j] = length_b2;
	}
	for (int i = 0; i < convex_a.size(); i++)
	{
		for (int j = 0; j < convex_b.size(); j++)
		{
			common_value[i][j] = -2 * (convex_a[i].x * convex_b[j].x + convex_a[i].y * convex_b[j].y);
		}
	}
}

double MovementArbitration::QueryTheMinkowskiSpace(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double theta, bool spin = false)
{
	/* 作最短距离排序 */
	std::priority_queue<MinkowsikiSpace, std::vector<MinkowsikiSpace>, CompareDistance> minimum_distance_to_origin;
	if (!spin) // 主旋最短距离计算
	{
		for (int i = 0; i < list_len_a.size(); i++)
		{
			for (int j = 0; j < list_len_b.size(); j++)
			{
				MinkowsikiSpace ms_st;
				double x = list_len_a[i] + list_len_b[j] + combine4cos[i][j] * cos(theta) + combine4sin[i][j] * sin(theta);
				ms_st.distance = sqrt(x);
				ms_st.i = i;
				ms_st.j = j;
				minimum_distance_to_origin.push(ms_st);
			}
		}
	}
	else // 自旋最短距离计算
	{
		for (int i = 0; i < spin_len_a.size(); i++)
		{
			for (int j = 0; j < spin_len_b.size(); j++)
			{
				MinkowsikiSpace ms_st;
				double x = spin_len_a[i] + spin_len_b[j] + common_value[i][j];
				ms_st.distance = sqrt(x);
				ms_st.i = i;
				ms_st.j = j;
				minimum_distance_to_origin.push(ms_st);
			}
		}
	}

	/* 找到离原点最近的一个点一定是闵可夫斯基差凸包上的点 */
	MinkowsikiSpace first = minimum_distance_to_origin.top();
	minimum_distance_to_origin.pop();

	/* 还原Minkowsiki空间向量 */
	Point first_a(convex_a[first.i].x * cos(theta) - convex_a[first.i].y * sin(theta), convex_a[first.i].y * cos(theta) + convex_a[first.i].x * sin(theta)); // 对convex_a点集进行旋转变换

	/* 最短的闵可夫斯基差向量一定在其凸包上 */
	Point vec_1(first_a.x - convex_b[first.j].x, first_a.y - convex_b[first.j].y);

	/* 遍历优先队列，不断更新vec_2使得获取与原点最近的距离*/

	double min_distance = first.distance; // 初始化最短距离为第一个点到原点的距离
	while (!minimum_distance_to_origin.empty())
	{
		MinkowsikiSpace second = minimum_distance_to_origin.top();
		minimum_distance_to_origin.pop();
		Point second_a(convex_a[second.i].x * cos(theta) - convex_a[second.i].y * sin(theta), convex_a[second.i].y * cos(theta) + convex_a[second.i].x * sin(theta));
		Point vec_2(second_a.x - convex_b[second.j].x, second_a.y - convex_b[second.j].y);
		Point seg(vec_2.x - vec_1.x, vec_2.y - vec_1.y);
		Point negative_seg(-1 * seg.x, -1 * seg.y);
		double length_seg = sqrt(seg.x * seg.x + seg.y * seg.y);
		double length_vec1_prj_seg = abs(ProjectOn(vec_1, seg) / length_seg);

		/*
		 * 用点乘判断原点是否在线段外取得最短距离
		 * 注意seg向量是vec_1指向vec_2的
		 */


		 /* 作vec_1与seg的点乘 */
		double prj_vec1_seg = ProjectOn(vec_1, negative_seg) / length_seg;
		/* 再作vec_2与negative_seg的点乘 */
		double prj_vec2_nseg = ProjectOn(vec_2, seg) / length_seg;

		/* 如果投影<0说明在线段外，如果都大于等于0，则用叉乘算平行四边形除以底边求高解决 */
		if (prj_vec1_seg < 0)
			min_distance = std::min(min_distance, first.distance);
		else if (prj_vec2_nseg < 0)
			min_distance = std::min(min_distance, second.distance);
		else
			min_distance = std::min(min_distance, abs(Cross(vec_1, vec_2) / length_seg));
	}
	return min_distance;
}

std::tuple<double, double> MovementArbitration::MainRotateApproach(const std::vector<Point>& detector, const std::vector<Point>& convex)
{
	Point origin(0, 0);
	Point min_vec = closest_point_to_origin(detector, origin);
	double min_len = sqrt(min_vec.x * min_vec.x + min_vec.y * min_vec.y);

	/* 预处理,作明可夫斯基差空间 */
	MinkowskiDifferenceSpace_Rotate(detector, convex);

	/* 用夹角求出最短顶点向量介于哪两个凸包点之间, 同时打表每个轮廓向量是否大于最短顶点距离*/
	double left_vec_cos_value = -2, right_vec_cos_value = -2;
	int left_vec_index = 0, right_vec_index = 0;
	for (int i = 0; i < convex.size(); i++)
	{
		/* 计算当前convex的长度 */
		double d = sqrt(convex[i].x * convex[i].x + convex[i].y * convex[i].y);

		/* 计算当前convex与min_vec的sin值，值越靠近0，向量间越靠近，正的在左边，负的在右边 */
		double sin_value = Cross(min_vec, convex[i]) / (min_len * d);
		double cos_value = ProjectOn(min_vec, convex[i]) / (min_len * d);
		if (sin_value > 0 && cos_value > left_vec_cos_value && d >= min_len)
		{
			left_vec_cos_value = cos_value;
			left_vec_index = i;
		}
		else if (sin_value < 0 && cos_value > right_vec_cos_value && d >= min_len)
		{
			right_vec_cos_value = cos_value;
			right_vec_index = i;
		}
		else if (sin_value == 0)
		{
			/* 如果等于0，检查一下是0°还是180°*/
			if (cos_value < 1e-16)
				left_vec_index = i; //如果0°更新左边界，因为默认是逆时针存储的
		}
	}

	/** 计算最左触礁边斜率，左边选x值较小的 凸包点逆时针存储，寻找上一个-1要注意加上容量再模除 **/
	double k2_left = (convex[(left_vec_index + convex.size() - 1) % convex.size()].y - convex[left_vec_index].y) / (convex[(left_vec_index + convex.size() - 1) % convex.size()].x - convex[left_vec_index].x);
	double x_left = -1 * (k2_left * convex[left_vec_index].y - k2_left * k2_left * convex[left_vec_index].x + sqrt(2 * convex[left_vec_index].x * convex[left_vec_index].y * k2_left - convex[left_vec_index].y * convex[left_vec_index].y - convex[left_vec_index].x * convex[left_vec_index].x * k2_left * k2_left + (1 + k2_left * k2_left) * min_len * min_len)) / (1 + k2_left * k2_left);
	double y_left = (convex[left_vec_index].y - k2_left * (convex[left_vec_index].x + sqrt(2 * convex[left_vec_index].x * convex[left_vec_index].y * k2_left - convex[left_vec_index].y * convex[left_vec_index].y - convex[left_vec_index].x * convex[left_vec_index].x * k2_left * k2_left + (1 + k2_left * k2_left) * min_len * min_len))) / (1 + k2_left * k2_left);


	/* 计算触礁向量与min_vec的夹角 */
	Point vec_left_strike(x_left, y_left);
	double length_left_strike = sqrt(vec_left_strike.x * vec_left_strike.x + vec_left_strike.y * vec_left_strike.y);
	double theta_left = acos(ProjectOn(vec_left_strike, min_vec) / (length_left_strike * min_len));

	/** 计算最右触礁边斜率，右边选x值较大的 **/
	double k2_right = (convex[(right_vec_index + 1) % convex.size()].y - convex[right_vec_index].y) / (convex[(right_vec_index + 1) % convex.size()].x - convex[right_vec_index].x);
	double x_right = -1 * (k2_right * convex[right_vec_index].y - k2_right * k2_right * convex[right_vec_index].x - sqrt(2 * convex[right_vec_index].x * convex[right_vec_index].y * k2_right - convex[right_vec_index].y * convex[right_vec_index].y - convex[right_vec_index].x * convex[right_vec_index].x * k2_right * k2_right + (1 + k2_right * k2_right) * min_len * min_len)) / (1 + k2_right * k2_right);
	double y_right = (convex[right_vec_index].y - k2_right * (convex[right_vec_index].x - sqrt(2 * convex[right_vec_index].x * convex[right_vec_index].y * k2_right - convex[right_vec_index].y * convex[right_vec_index].y - convex[right_vec_index].x * convex[right_vec_index].x * k2_right * k2_right + (1 + k2_right * k2_right) * min_len * min_len))) / (1 + k2_right * k2_right);


	/* 计算左触礁向量与min_vec的夹角 */
	Point vec_right_strike(x_right, y_right);
	double length_right_strike = sqrt(vec_right_strike.x * vec_right_strike.x + vec_right_strike.y * vec_right_strike.y);
	double theta_right = acos(ProjectOn(vec_right_strike, min_vec) / (length_right_strike * min_len));


	/* 步进查询，步进值为1度 */
	double left_max = FindMaxCCWAngle(detector, convex, theta_left, M_PI * 0.01 / 180);
	double right_max = FindMaxCWAngle(detector, convex, theta_right, M_PI * 0.01 / 180);
	return { left_max, right_max };
}

double MovementArbitration::ProjectOn(Point vec1, Point vec2)
{
	return vec1.x * vec2.x + vec1.y * vec2.y;
}

double MovementArbitration::Cross(Point vec1, Point vec2)
{
	return vec1.x * vec2.y - vec1.y * vec2.x;
}

double MovementArbitration::dist(Point p1, Point p2)
{
	return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

Point MovementArbitration::closest_point_on_segment(Point p1, Point p2, Point origin)
{
	double l2 = pow(dist(p1, p2), 2);
	if (l2 == 0.0) return p1;
	double t = std::max(0.0, std::min(1.0, ((origin.x - p1.x) * (p2.x - p1.x) + (origin.y - p1.y) * (p2.y - p1.y)) / l2));
	Point projection = Point(p1.x + t * (p2.x - p1.x), p1.y + t * (p2.y - p1.y));
	return projection;
}

Point MovementArbitration::closest_point_to_origin(std::vector<Point> point_list, Point origin)
{
	double min_dist = dist(point_list[0], origin);
	Point closest_point = point_list[0];
	for (size_t i = 0; i < point_list.size(); ++i) {
		Point p1 = point_list[i];
		Point p2 = point_list[(i + 1) % point_list.size()];
		Point new_point = closest_point_on_segment(p1, p2, origin);
		double new_dist = dist(new_point, origin);
		if (new_dist < min_dist) {
			min_dist = new_dist;
			closest_point = new_point;
		}
	}
	return closest_point;
}

Point MovementArbitration::RectangleCenter(const std::vector<Point>& rectangle)
{
	double sum_x = 0, sum_y = 0;
	for (const Point& p : rectangle) {
		sum_x += p.x;
		sum_y += p.y;
	}
	return { sum_x / 4, sum_y / 4 };
}

Point MovementArbitration::RotateCCW(const Point& p, const Point& center, double theta)
{
	double x_new = (p.x - center.x) * cos(theta) - (p.y - center.y) * sin(theta) + center.x;
	double y_new = (p.x - center.x) * sin(theta) + (p.y - center.y) * cos(theta) + center.y;
	return { x_new, y_new };
}

std::vector<Point> MovementArbitration::RotateRectangle(const std::vector<Point>& rectangle, double theta)
{
	Point center = RectangleCenter(rectangle);
	std::vector<Point> rotated_rectangle;
	for (const Point& p : rectangle) {
		rotated_rectangle.push_back(RotateCCW(p, center, theta));
	}
	return rotated_rectangle;
}

double MovementArbitration::FindMaxCCWAngle(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double step)
{
	double max_angle = 0, r = 0;
	bool found = false, dmodel = false;
	MinkowskiDifferenceSpace_Spin(convex_a, convex_b, 0);
	r = QueryTheMinkowskiSpace(convex_a, convex_b, 0, true);
	if (r > 1.0)dmodel = true; // 如果大于1，则以降到1为目标
	else dmodel = false; // 否则以升到1为目标
	for (double angle = 0; angle <= 2 * M_PI; angle += step)
	{

		/* 本质上就是查询该探头状态下主旋转为0的距离情况 */
		auto rotated_detector = RotateRectangle(convex_a, angle);
		MinkowskiDifferenceSpace_Spin(rotated_detector, convex_b, angle);
		r = QueryTheMinkowskiSpace(rotated_detector, convex_b, 0, true);

		std::cout << "Left r: " << r << " L radian: " << angle << std::endl;
		if (dmodel)
		{
			if (r - 1.0 > EPS)
			{
				max_angle = angle;
			}
			else
			{
				found = true;
				break;
			}
		}
		else
		{
			if (1.0 - r > EPS)
			{
				max_angle = angle;
			}
			else
			{
				found = true;
				break;
			}
		}
	}
	return found ? max_angle : 0;
}

double MovementArbitration::FindMaxCWAngle(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double step)
{
	double max_angle = 0, r = 0;
	bool found = false, dmodel = false;
	MinkowskiDifferenceSpace_Spin(convex_a, convex_b, 0);
	r = QueryTheMinkowskiSpace(convex_a, convex_b, 0, true);
	if (r > 1.0)dmodel = true; // 如果大于1，则以降到1为目标
	else dmodel = false; // 否则以升到1为目标
	for (double angle = 0; angle <= 2 * M_PI; angle += step)
	{

		/* 本质上就是查询该探头状态下主旋转为0的距离情况 */
		auto rotated_detector = RotateRectangle(convex_a, -angle);
		MinkowskiDifferenceSpace_Spin(rotated_detector, convex_b, -angle);
		r = QueryTheMinkowskiSpace(rotated_detector, convex_b, 0, true);
		std::cout << "Right r: " << r << " R radian: " << angle << std::endl;
		if (dmodel)
		{
			if (r - 1.0 > EPS)
			{
				max_angle = angle;
			}
			else
			{
				found = true;
				break;
			}
		}
		else
		{
			if (1.0 - r > EPS)
			{
				max_angle = angle;
			}
			else
			{
				found = true;
				break;
			}
		}
	}
	return found ? max_angle : 0;
}

double MovementArbitration::FindMaxCWAngle(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double angle_right, double step)
{
	double max_angle = 0, r = 0;
	bool found = false, dmodel = false;
	r = QueryTheMinkowskiSpace(convex_a, convex_b, 0);
	if (r > 1.0)dmodel = true; // 如果大于1，则以降到1为目标
	else dmodel = false; // 否则以升到1为目标
	for (double angle = 0; angle >= -angle_right; angle -= step)
	{
		r = QueryTheMinkowskiSpace(convex_a, convex_b, angle);
		//std::cout << "Right r: " << r << " R radian: " << angle << std::endl;
		if (dmodel)
		{
			if (r - 1.0 > EPS)
			{
				max_angle = angle;
			}
			else
			{
				found = true;
				break;
			}
		}
		else
		{
			if (1.0 - r > EPS)
			{
				max_angle = angle;
			}
			else
			{
				found = true;
				break;
			}
		}
	}
	return found ? max_angle : 0;
}

double MovementArbitration::FindMaxCCWAngle(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double angle_left, double step)
{
	double max_angle = 0, r = 0;
	bool found = false, dmodel = false;
	r = QueryTheMinkowskiSpace(convex_a, convex_b, 0);
	if (r > 1.0)dmodel = true; // 如果大于1，则以降到1为目标
	else dmodel = false; // 否则以升到1为目标
	for (double angle = 0; angle <= angle_left; angle += step)
	{
		r = QueryTheMinkowskiSpace(convex_a, convex_b, angle);
		//std::cout << "Left r: " << r << " L radian: " << angle << std::endl;
		if (dmodel)
		{
			if (r - 1.0 > EPS)
			{
				max_angle = angle;
			}
			else
			{
				found = true;
				break;
			}
		}
		else
		{
			if (1.0 - r > EPS)
			{
				max_angle = angle;
			}
			else
			{
				found = true;
				break;
			}
		}
	}
	return found ? max_angle : 0;
}