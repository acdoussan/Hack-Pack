// Adam Doussan AD844156 04/14/2017

public class LinePlaneIntersect
{
	// example where input is given as points
	public static void main(String [] args)
	{
		Scanner in = new Scanner(System.in);

		int run = in.nextInt();

		for(int r = 1; r <= run; r++)
		{
			pt p1 = new pt(in.nextInt(),in.nextInt(),in.nextInt());
			pt p2 = new pt(in.nextInt(),in.nextInt(),in.nextInt());

			line l = new line(p1.vect(p2), p1);

			pt p3 = new pt(in.nextInt(),in.nextInt(),in.nextInt());
			pt p4 = new pt(in.nextInt(),in.nextInt(),in.nextInt());
			pt p5 = new pt(in.nextInt(),in.nextInt(),in.nextInt());

			plane p = new plane(p3.vect(p4), p3.vect(p5), p3);

			System.out.format("Data Set #%d:\n", r);
			solve(l,p);
			System.out.println();
		}
	}

	public static void solve(line l, plane p)
	{
		int right = p.d - (l.p.x*p.orth.x + l.p.y*p.orth.y + l.p.z*p.orth.z);
		int left = (l.vect.x*p.orth.x + l.vect.y*p.orth.y + l.vect.z*p.orth.z);

		if(right == 0)
		{
			System.out.println("The line lies on the plane.");
		}

		else if(left == 0)
		{
			System.out.println("There is no intersection.");
		}

		else
		{
			double lam = ((double)right) / left;
			System.out.format("The intersection is the point (%.1f, %.1f, %.1f).\n",
			l.p.x + l.vect.x*lam, l.p.y + l.vect.y*lam, l.p.z + l.vect.z*lam);
		}
	}
}

class pt
{
	int x, y, z;

	pt(int x, int y, int z)
	{
		this.x = x; this.y = y; this.z = z;
	}

	public pt vect(pt pt2)
	{
		return new pt(pt2.x-x, pt2.y-y, pt2.z-z);
	}

	public pt cross(pt pt2)
	{
		int x = this.y * pt2.z - pt2.y * this.z;
		int y = this.z * pt2.x - pt2.z * this.x;
		int z = this.x * pt2.y - pt2.x * this.y;

		return new pt(x,y,z);
	}

	public int solve(pt pt2)
	{
		return this.x*pt2.x + this.y*pt2.y + this.z*pt2.z;
	}
}

class line
{
	pt vect, p;

	public line(pt vect, pt p)
	{
		this.vect = vect; this.p = p;
	}
}

class plane
{
	pt orth;
	int d;

	public plane(pt vect1, pt vect2, pt ori)
	{
		orth = vect1.cross(vect2);
		d = orth.solve(ori);
	}
}

