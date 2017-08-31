// Adam Doussan 4/24/2017

public class LineLineIntersection
{
	// java.awt.geom.Line2D.linesIntersect()

	public static boolean intersect(line a, line b)
	{
		return a.intersect(b);
	}
}

class vect
{
	public double x,y;

	public vect(pt start, pt end)
	{
		x = end.x - start.x;
		y = end.y - start.y;
	}
}

class pt
{
	public double x, y;

	public pt(double x, double y)
	{
		this.x = x; this.y = y;
	}
}

class line
{
	public pt p;
	public vect dir;

	public line(pt start, pt end)
	{
		p = start;
		dir = new vect(start, end);
	}

	public boolean intersect(line other)
	{
		double den = det(dir.x, -other.dir.x, dir.y, -other.dir.y);
		double numLam = det(other.p.x-p.x, -other.dir.x, other.p.y-p.y, -other.dir.y);

		// paralell
		if(den < 1e-9)
			return false;

		else if(numLam/den > 1)
			return false;
		
		else
			return true;
	}

	public static double det(double a,double b ,double c ,double d)
	{
		return a*d-b*c;
	}
}
