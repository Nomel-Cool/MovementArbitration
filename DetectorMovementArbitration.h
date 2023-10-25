#pragma once
#ifndef DETECT_MOVEMENT_ARBITRATION
#define DETECT_MOVEMENT_ARBITRATION
#include <math.h>
#include <tuple>
#include <vector>
#include <queue>

#include <iostream>

constexpr double M_PI = 3.14159265358979323846;
constexpr double EPS = 1e-2; // 迫近精度
/*TODO（有生之年）
 代码整理工作：
	1. 把自旋和主转的闵可夫斯基差空间计算分开
	2. 把自旋和主旋的查询闵可夫斯基差空间分开
	3. 优化自旋的闵可夫斯基差空间计算，直接加入旋转后做公式常量分解计算，这样可以和主旋的查询统一格式
*/


/// <summary>
/// 可以证明：闵可夫斯基差空间里面模长最短的两个向量，一定就是整个闵可夫斯基凸包离原点最短的线段
/// 最短的点一定是凸包上的点，倒数第二短的点一定也在凸包上且与最短点共线，由凸包性质反证可知如果不共线，一定存在比它更短的点与最短点共线，这与设定矛盾
/// </summary>
struct MinkowsikiSpace
{
	double distance;
	int i, j; //i和j分别对应两个凸包的两个点
};

struct CompareDistance {
	bool operator()(MinkowsikiSpace const& p1, MinkowsikiSpace const& p2) {
		// 返回true的元素排在后面
		return p1.distance > p2.distance;
	}
};

struct Point
{
	double x, y;
	Point() :x(0), y(0) {}
	Point(double _x, double _y) { x = _x; y = _y; }
};
/// <summary>
/// 任务 for D:
/// 输入：
/// 1. 防碰撞算法给出的最短碰撞距离d
/// 2. 探头规模
/// 3. 探头中心点坐标
/// 4. 探头自转角度
/// 5. 探头主转角度
/// 
/// 目标：
/// 1. 探头离人体最近1CM
/// 
/// 处理：
/// 1. 根据当前探头状态和可到达的运动状态判断返回是否可以使探头离近
/// 2. 给出离近的运动探头的控制方程ξ(φ,θ,ρ)，其中φ是探头自转角度，θ是主转角度，ρ是控制距离，它是ρ=|d - 1|，d是探头当前离人距离，可正，可负
///    当d = 1即ρ = 0时，φ和θ不需要改变即为解
///    不妨设F(φ,θ)=ξ(φ,θ,0)，即保证当ρ=0时，φ和θ保持控制恒等的所有解的范围
///    这意味着需要找到一个初始解使得ρ=0
/// 3. 由于对象被转化为凸包，对凸包外接一个1CM距离的凸包多边形，那么整个外接凸包多边形都是解，从机械角度来看，可以认为d始终大于1
/// 4. 那么探头可以想象成一条逼近安全凸多边形 C 的直线 L
/// 5. 直线 L接触安全凸多边形 C 的方式有：
/// 5.1 d维度直接逼近触碰，第一种可能：先会碰到点，第二种可能：先碰到边
/// 5.2 自旋使得边界触碰
/// 5.3 如果主旋中心在轮廓内，求轮廓向量最短值 a，以及探头平面的径向距离 b，若 a ≥ b，则一定可以使探头主旋使触碰，若a < b，则一定不可能触碰（因为探头初始值就是安全的）
/// 5.4 如果主旋中心在轮廓外，也是同样的判断
/// 
/// 控制论上可以对这三个操作做偏导，也就是这三者可以独立判断
/// </summary>
class MovementArbitration
{
public:
	/// <summary>
	/// 径向到达
	/// </summary>
	/// <param name="vec_detector">探头向量</param>
	/// <param name="approach_point">可能是探头，也可能是人体轮廓的接近点</param>
	/// <param name="start_point">可能是探头，也可能是人体轮廓的接近边向量的起始点</param>
	/// <param name="end_point">可能是探头，也可能是人体轮廓的接近边向量的终点</param>
	/// <param name="safe_distance">项目要求安全距离</param>
	/// <returns>径向移动后的探头位置</returns>
	Point RadiusDirectionApproach(const Point& vec_detector, const Point& approach_point, const Point& start_point, const Point& end_point, double safe_distance = 1);

	/// <summary>
	/// 自旋到达
	/// </summary>
	/// <param name="detector">探头向量</param>
	/// <param name="convex">轮廓向量集</param>
	/// <returns>顺时针，逆时针旋转角度容差</returns>
	std::tuple<double, double> Spin2Approach(const std::vector<Point>& detector, const std::vector<Point>& convex);

	/// <summary>
	/// 主旋转明可夫斯基差空间，利用优先队列nlogn的排序时间复杂度，对化简后m×n个点的距离进行排序，并挑出最近的两个点作线段求原点到该线段的距离
	/// 注意这里的打表是用矩形的点减去轮廓的点构建的闵可夫斯基差向量
	/// </summary>
	/// <param name="convex_a">矩形凸包点集，转</param>
	/// <param name="convex_b">轮廓凸包点集,不转</param>
	void MinkowskiDifferenceSpace_Rotate(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b);

	/// <summary>
	/// 自转明可夫斯基差空间，利用优先队列nlogn的排序时间复杂度，对化简后m×n个点的距离进行排序，并挑出最近的两个点作线段求原点到该线段的距离
	/// 注意这里的打表是用矩形的点减去轮廓的点构建的闵可夫斯基差向量
	/// </summary>
	/// <param name="convex_a">矩形凸包点集，转</param>
	/// <param name="convex_b">轮廓凸包点集,不转</param>
	/// <param name="theta">矩形凸包点集自转角度</param>
	void MinkowskiDifferenceSpace_Spin(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double theta);

	/// <summary>
	/// 查询逆时针旋转theta度后两个凸包的距离
	/// </summary>
	/// <param name="convex_a">矩形凸包点集，转</param>
	/// <param name="convex_b">轮廓凸包点集,不转</param>
	/// <param name="theta">逆时针旋转角度</param>
	/// <returns>两个凸多边形最短距离</returns>
	double QueryTheMinkowskiSpace(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double theta, bool spin = false);

	/// <summary>
	/// 主旋到达
	/// </summary>
	/// <param name="detector">探头向量</param>
	/// <param name="convex">轮廓向量集</param>
	/// <returns>顺时针，逆时针旋转角度容差</returns>
	std::tuple<double, double> MainRotateApproach(const std::vector<Point>& detector, const std::vector<Point>& convex);


private:

	std::vector<double> list_len_a, // 记录 ux^2 + uy^2
						list_len_b, // 记录 vx^2 + vy^2
						spin_len_a, // 记录 自旋后的 ux^2 + uy^2
						spin_len_b; // 记录 自旋后的 vx^2 + vy^2

	std::vector<std::vector<double>> combine4cos, // 记录主旋明可夫斯基距离cos的系数
									 combine4sin, // 记录主旋明可夫斯基距离sin的系数
									 common_value; // 记录自旋旋明可夫斯基距离的常量


	/// <summary>
	/// vec1到vec2的投影值
	/// </summary>
	/// <param name="vec1">投影向量</param>
	/// <param name="vec2">背投影向量</param>
	/// <returns>投影值</returns>
	double ProjectOn(Point vec1, Point vec2);

	/// <summary>
	/// 叉乘
	/// </summary>
	/// <param name="vec1"></param>
	/// <param name="vec2"></param>
	/// <returns></returns>
	double Cross(Point vec1, Point vec2);

	/// <summary>
	/// 点距
	/// </summary>
	/// <param name="p1"></param>
	/// <param name="p2"></param>
	/// <returns></returns>
	double dist(Point p1, Point p2);

	/// <summary>
	/// 点origin到segment(p1,p2)的最近距离
	/// </summary>
	/// <param name="p1"></param>
	/// <param name="p2"></param>
	/// <param name="origin"></param>
	/// <returns></returns>
	Point closest_point_on_segment(Point p1, Point p2, Point origin);

	/// <summary>
	/// 点集到原点的最近距离
	/// </summary>
	/// <param name="point_list"></param>
	/// <param name="origin"></param>
	/// <returns></returns>
	Point closest_point_to_origin(std::vector<Point> point_list, Point origin);

	/// <summary>
	/// 获取给定矩形的形心坐标
	/// </summary>
	/// <param name="rectangle"></param>
	/// <returns></returns>
	Point RectangleCenter(const std::vector<Point>& rectangle);

	/// <summary>
	/// 逆时针旋转theta弧度
	/// </summary>
	/// <param name="p"></param>
	/// <param name="center"></param>
	/// <param name="theta"></param>
	/// <returns></returns>
	Point RotateCCW(const Point& p, const Point& center, double theta);

	/// <summary>
	/// 对给定矩形逆时针旋转theta弧度
	/// </summary>
	/// <param name="rectangle"></param>
	/// <param name="theta"></param>
	/// <returns></returns>
	std::vector<Point> RotateRectangle(const std::vector<Point>& rectangle, double theta);

	/// <summary>
	/// (重载)线性寻找自旋过程中离安全线最近的逆旋转角度
	/// </summary>
	/// <param name="convex_a"></param>
	/// <param name="convex_b"></param>
	/// <param name="step"></param>
	/// <returns></returns>
	double FindMaxCCWAngle(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double step);

	/// <summary>
	/// (重载)线性寻找自旋过程中离安全线最近的顺旋转角度
	/// </summary>
	/// <param name="convex_a"></param>
	/// <param name="convex_b"></param>
	/// <param name="step"></param>
	/// <returns></returns>
	double FindMaxCWAngle(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double step);

	/// <summary>
	/// 线性寻找主旋过程中离安全线最近的顺旋转角度
	/// </summary>
	/// <param name="convex_a"></param>
	/// <param name="convex_b"></param>
	/// <param name="angle_right"></param>
	/// <param name="step"></param>
	/// <returns></returns>
	double FindMaxCWAngle(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double angle_right, double step);

	/// <summary>
	/// 线性寻找主旋过程中离安全线最近的逆旋转角度
	/// </summary>
	/// <param name="convex_a"></param>
	/// <param name="convex_b"></param>
	/// <param name="angle_left"></param>
	/// <param name="step"></param>
	/// <returns></returns>
	double FindMaxCCWAngle(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double angle_left, double step);
};

#endif // !DETECT_MOVEMENT_ARBITRATION
