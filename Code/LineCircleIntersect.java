// Ross 4/24/2017

public class LineCircleIntersect
{
	// returns true if the segment intersects the circle, false otherwise
    public static boolean intersect(pt start, pt end, circ circle)
    {
		double dx = end.x - start.x;
		double dy = end.y - start.y;

		double a = ((dx*dx)+(dy*dy));
		double b = 2 * (dx * (start.x - circle.x) + dy * (start.y - circle.y));
		double c = ((start.x - circle.x) * (start.x - circle.x))
					+ ((start.y - circle.y) * (start.y - circle.y))
					- (circle.radius * circle.radius);

		double det = (b*b) - (4 * a * c);

		// Answer will be non-real; discard.
		if (det < 0)
			return false;

		double soln1 = ((-1 * b) + Math.sqrt(det))/(2*a);
		double soln2 = ((-1 * b) - Math.sqrt(det))/(2*a);

		if (((soln1 < 0) || (soln1 > 1)) && ((soln2 < 0) || (soln2 > 1)))
			return false;
		return true;
    }
}

class pt
{
	double x,y;
}

class circ
{
	double x,y;
	double radius;
}
