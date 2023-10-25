#pragma once
#ifndef DETECT_MOVEMENT_ARBITRATION
#define DETECT_MOVEMENT_ARBITRATION
#include <math.h>
#include <tuple>
#include <vector>
#include <queue>

#include <iostream>

constexpr double M_PI = 3.14159265358979323846;
constexpr double EPS = 1e-2; // �Ƚ�����
/*TODO������֮�꣩
 ������������
	1. ����������ת���ɿɷ�˹����ռ����ֿ�
	2. �������������Ĳ�ѯ�ɿɷ�˹����ռ�ֿ�
	3. �Ż��������ɿɷ�˹����ռ���㣬ֱ�Ӽ�����ת������ʽ�����ֽ���㣬�������Ժ������Ĳ�ѯͳһ��ʽ
*/


/// <summary>
/// ����֤�����ɿɷ�˹����ռ�����ģ����̵�����������һ�����������ɿɷ�˹��͹����ԭ����̵��߶�
/// ��̵ĵ�һ����͹���ϵĵ㣬�����ڶ��̵ĵ�һ��Ҳ��͹����������̵㹲�ߣ���͹�����ʷ�֤��֪��������ߣ�һ�����ڱ������̵ĵ�����̵㹲�ߣ������趨ì��
/// </summary>
struct MinkowsikiSpace
{
	double distance;
	int i, j; //i��j�ֱ��Ӧ����͹����������
};

struct CompareDistance {
	bool operator()(MinkowsikiSpace const& p1, MinkowsikiSpace const& p2) {
		// ����true��Ԫ�����ں���
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
/// ���� for D:
/// ���룺
/// 1. ����ײ�㷨�����������ײ����d
/// 2. ̽ͷ��ģ
/// 3. ̽ͷ���ĵ�����
/// 4. ̽ͷ��ת�Ƕ�
/// 5. ̽ͷ��ת�Ƕ�
/// 
/// Ŀ�꣺
/// 1. ̽ͷ���������1CM
/// 
/// ����
/// 1. ���ݵ�ǰ̽ͷ״̬�Ϳɵ�����˶�״̬�жϷ����Ƿ����ʹ̽ͷ���
/// 2. ����������˶�̽ͷ�Ŀ��Ʒ��̦�(��,��,��)�����Ц���̽ͷ��ת�Ƕȣ�������ת�Ƕȣ����ǿ��ƾ��룬���Ǧ�=|d - 1|��d��̽ͷ��ǰ���˾��룬�������ɸ�
///    ��d = 1���� = 0ʱ���պͦȲ���Ҫ�ı伴Ϊ��
///    ������F(��,��)=��(��,��,0)������֤����=0ʱ���պͦȱ��ֿ��ƺ�ȵ����н�ķ�Χ
///    ����ζ����Ҫ�ҵ�һ����ʼ��ʹ�æ�=0
/// 3. ���ڶ���ת��Ϊ͹������͹�����һ��1CM�����͹������Σ���ô�������͹������ζ��ǽ⣬�ӻ�е�Ƕ�������������Ϊdʼ�մ���1
/// 4. ��ô̽ͷ���������һ���ƽ���ȫ͹����� C ��ֱ�� L
/// 5. ֱ�� L�Ӵ���ȫ͹����� C �ķ�ʽ�У�
/// 5.1 dά��ֱ�ӱƽ���������һ�ֿ��ܣ��Ȼ������㣬�ڶ��ֿ��ܣ���������
/// 5.2 ����ʹ�ñ߽紥��
/// 5.3 ������������������ڣ��������������ֵ a���Լ�̽ͷƽ��ľ������ b���� a �� b����һ������ʹ̽ͷ����ʹ��������a < b����һ�������ܴ�������Ϊ̽ͷ��ʼֵ���ǰ�ȫ�ģ�
/// 5.4 ������������������⣬Ҳ��ͬ�����ж�
/// 
/// �������Ͽ��Զ�������������ƫ����Ҳ���������߿��Զ����ж�
/// </summary>
class MovementArbitration
{
public:
	/// <summary>
	/// ���򵽴�
	/// </summary>
	/// <param name="vec_detector">̽ͷ����</param>
	/// <param name="approach_point">������̽ͷ��Ҳ���������������Ľӽ���</param>
	/// <param name="start_point">������̽ͷ��Ҳ���������������Ľӽ�����������ʼ��</param>
	/// <param name="end_point">������̽ͷ��Ҳ���������������Ľӽ����������յ�</param>
	/// <param name="safe_distance">��ĿҪ��ȫ����</param>
	/// <returns>�����ƶ����̽ͷλ��</returns>
	Point RadiusDirectionApproach(const Point& vec_detector, const Point& approach_point, const Point& start_point, const Point& end_point, double safe_distance = 1);

	/// <summary>
	/// ��������
	/// </summary>
	/// <param name="detector">̽ͷ����</param>
	/// <param name="convex">����������</param>
	/// <returns>˳ʱ�룬��ʱ����ת�Ƕ��ݲ�</returns>
	std::tuple<double, double> Spin2Approach(const std::vector<Point>& detector, const std::vector<Point>& convex);

	/// <summary>
	/// ����ת���ɷ�˹����ռ䣬�������ȶ���nlogn������ʱ�临�Ӷȣ��Ի����m��n����ľ���������򣬲�������������������߶���ԭ�㵽���߶εľ���
	/// ע������Ĵ�����þ��εĵ��ȥ�����ĵ㹹�����ɿɷ�˹��������
	/// </summary>
	/// <param name="convex_a">����͹���㼯��ת</param>
	/// <param name="convex_b">����͹���㼯,��ת</param>
	void MinkowskiDifferenceSpace_Rotate(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b);

	/// <summary>
	/// ��ת���ɷ�˹����ռ䣬�������ȶ���nlogn������ʱ�临�Ӷȣ��Ի����m��n����ľ���������򣬲�������������������߶���ԭ�㵽���߶εľ���
	/// ע������Ĵ�����þ��εĵ��ȥ�����ĵ㹹�����ɿɷ�˹��������
	/// </summary>
	/// <param name="convex_a">����͹���㼯��ת</param>
	/// <param name="convex_b">����͹���㼯,��ת</param>
	/// <param name="theta">����͹���㼯��ת�Ƕ�</param>
	void MinkowskiDifferenceSpace_Spin(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double theta);

	/// <summary>
	/// ��ѯ��ʱ����תtheta�Ⱥ�����͹���ľ���
	/// </summary>
	/// <param name="convex_a">����͹���㼯��ת</param>
	/// <param name="convex_b">����͹���㼯,��ת</param>
	/// <param name="theta">��ʱ����ת�Ƕ�</param>
	/// <returns>����͹�������̾���</returns>
	double QueryTheMinkowskiSpace(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double theta, bool spin = false);

	/// <summary>
	/// ��������
	/// </summary>
	/// <param name="detector">̽ͷ����</param>
	/// <param name="convex">����������</param>
	/// <returns>˳ʱ�룬��ʱ����ת�Ƕ��ݲ�</returns>
	std::tuple<double, double> MainRotateApproach(const std::vector<Point>& detector, const std::vector<Point>& convex);


private:

	std::vector<double> list_len_a, // ��¼ ux^2 + uy^2
						list_len_b, // ��¼ vx^2 + vy^2
						spin_len_a, // ��¼ ������� ux^2 + uy^2
						spin_len_b; // ��¼ ������� vx^2 + vy^2

	std::vector<std::vector<double>> combine4cos, // ��¼�������ɷ�˹������cos��ϵ��
									 combine4sin, // ��¼�������ɷ�˹������sin��ϵ��
									 common_value; // ��¼���������ɷ�˹������ĳ���


	/// <summary>
	/// vec1��vec2��ͶӰֵ
	/// </summary>
	/// <param name="vec1">ͶӰ����</param>
	/// <param name="vec2">��ͶӰ����</param>
	/// <returns>ͶӰֵ</returns>
	double ProjectOn(Point vec1, Point vec2);

	/// <summary>
	/// ���
	/// </summary>
	/// <param name="vec1"></param>
	/// <param name="vec2"></param>
	/// <returns></returns>
	double Cross(Point vec1, Point vec2);

	/// <summary>
	/// ���
	/// </summary>
	/// <param name="p1"></param>
	/// <param name="p2"></param>
	/// <returns></returns>
	double dist(Point p1, Point p2);

	/// <summary>
	/// ��origin��segment(p1,p2)���������
	/// </summary>
	/// <param name="p1"></param>
	/// <param name="p2"></param>
	/// <param name="origin"></param>
	/// <returns></returns>
	Point closest_point_on_segment(Point p1, Point p2, Point origin);

	/// <summary>
	/// �㼯��ԭ����������
	/// </summary>
	/// <param name="point_list"></param>
	/// <param name="origin"></param>
	/// <returns></returns>
	Point closest_point_to_origin(std::vector<Point> point_list, Point origin);

	/// <summary>
	/// ��ȡ�������ε���������
	/// </summary>
	/// <param name="rectangle"></param>
	/// <returns></returns>
	Point RectangleCenter(const std::vector<Point>& rectangle);

	/// <summary>
	/// ��ʱ����תtheta����
	/// </summary>
	/// <param name="p"></param>
	/// <param name="center"></param>
	/// <param name="theta"></param>
	/// <returns></returns>
	Point RotateCCW(const Point& p, const Point& center, double theta);

	/// <summary>
	/// �Ը���������ʱ����תtheta����
	/// </summary>
	/// <param name="rectangle"></param>
	/// <param name="theta"></param>
	/// <returns></returns>
	std::vector<Point> RotateRectangle(const std::vector<Point>& rectangle, double theta);

	/// <summary>
	/// (����)����Ѱ�������������밲ȫ�����������ת�Ƕ�
	/// </summary>
	/// <param name="convex_a"></param>
	/// <param name="convex_b"></param>
	/// <param name="step"></param>
	/// <returns></returns>
	double FindMaxCCWAngle(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double step);

	/// <summary>
	/// (����)����Ѱ�������������밲ȫ�������˳��ת�Ƕ�
	/// </summary>
	/// <param name="convex_a"></param>
	/// <param name="convex_b"></param>
	/// <param name="step"></param>
	/// <returns></returns>
	double FindMaxCWAngle(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double step);

	/// <summary>
	/// ����Ѱ�������������밲ȫ�������˳��ת�Ƕ�
	/// </summary>
	/// <param name="convex_a"></param>
	/// <param name="convex_b"></param>
	/// <param name="angle_right"></param>
	/// <param name="step"></param>
	/// <returns></returns>
	double FindMaxCWAngle(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double angle_right, double step);

	/// <summary>
	/// ����Ѱ�������������밲ȫ�����������ת�Ƕ�
	/// </summary>
	/// <param name="convex_a"></param>
	/// <param name="convex_b"></param>
	/// <param name="angle_left"></param>
	/// <param name="step"></param>
	/// <returns></returns>
	double FindMaxCCWAngle(const std::vector<Point>& convex_a, const std::vector<Point>& convex_b, double angle_left, double step);
};

#endif // !DETECT_MOVEMENT_ARBITRATION
