// Adam Doussan AD844156 04/28/2017

public class Polygon
{
	public static boolean ptInPoly(pt intersect, pt [] poly)
	{
		double totalAngle = 0;

		for(int i = 0; i < poly.length; i++)
		{
			pt v1 = intersect.vect(poly[i]);
			pt v2 = intersect.vect(poly[(i+1)%poly.length]);

			try
			{
				totalAngle += Math.acos((v1.dot(v2)) / (v1.mag() * v2.mag()));
			}
			catch(Exception e)
			{
			}
		}

		return (2*Math.PI) - totalAngle < 1e-9;
	}

	// assumes polygon is 2d, z is unused and set to 0, and points are ordered
	// such that pt[x] -> pt[(x+1) % pt.length] is connected with a line. Also
	// assumes simple polygon
	public static double area(pt [] poly)
	{
		double ans = 0;
		pt ori = new pt(0,0,0);

		for(int i = 0; i < poly.length; i++)
		{
			pt v1 = ori.vect(poly[i]);
			pt v2 = ori.vect(poly[(i+1)%poly.length]);
			ans += v1.cross(v2).z;
		}

		return 0.5 * Math.abs(ans);
	}
}

class pt
{
	public double x, y, z;

	pt(double x, double y, double z)
	{
		this.x = x; this.y = y; this.z = z;
	}

	public pt vect(pt pt2)
	{
		return new pt(pt2.x-x, pt2.y-y, pt2.z-z);
	}

	public pt cross(pt pt2)
	{
		double x = this.y * pt2.z - pt2.y * this.z;
		double y = this.z * pt2.x - pt2.z * this.x;
		double z = this.x * pt2.y - pt2.x * this.y;

		return new pt(x,y,z);
	}

	public double mag()
	{
		return Math.sqrt(Math.pow(this.x,2) + Math.pow(this.y,2) + Math.pow(this.z,2));
	}

	public double dot(pt vect)
	{
		return (this.x * vect.x) + (this.y * vect.y) + (this.z * vect.z);
	}
}
