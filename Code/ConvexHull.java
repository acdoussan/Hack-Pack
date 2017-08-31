// Adam Doussan AD844156 04/15/2017

public class ConvexHull
{
	public static void main(String [] args)
	{
		Scanner in = new Scanner(System.in);
		int cases = in.nextInt();

		for(int i = 0; i < cases; i++)
		{
			int numCase = in.nextInt();
			int numPoints = in.nextInt();
			pt[] pts = new pt[numPoints];

			for(int k = 0; k < numPoints; k++)
				pts[k] = new pt(in.nextInt(), in.nextInt());

			pt.ref = getRef(pts, numPoints);

			ArrayList<pt> ans = grahamScan(pts, numPoints);
		}
	}

	// returns the point in pts with min y breaking tie by min x
	public static pt getRef(pt[] pts, int n)
	{
		pt ans = pts[0];

		for(int i = 1; i < n; i++)	
			if(pts[i].y < ans.y || (pts[i].y == ans.y && pts[i].x < ans.x))			      
				ans = pts[i];

		return ans;
	}

	public static ArrayList<pt> grahamScan(pt[] pts, int n)
	{
		Arrays.sort(pts);

		Stack<pt> st = new Stack<>();
		st.push(pts[0]);
		st.push(pts[1]);	
	
		for(int i = 2; i < n; i++)
		{
			pt next = pts[i];
			pt mid = st.pop();
			pt prev = st.pop();

			while(!prev.isRightTurn(mid, next))
			{
				mid = prev;
				prev = st.pop();
			}

			st.push(prev);
			if(!prev.isStraight(mid, next)) // deleting collinear points
				st.push(mid);
			st.push(next);
		}
	
		ArrayList<pt> ans = new ArrayList<>();

		while(!st.isEmpty())
			ans.add(st.pop());

		return ans;
	}

}

class pt implements Comparable<pt>
{
	public static pt ref;
	public int x, y;

	public pt(int x, int y)
	{
		this.x = x; this.y = y;
	}

	public pt vect(pt pt2)
	{
		return new pt(pt2.x-x, pt2.y-y);
	}

	public double dist(pt pt2)
	{
		return Math.sqrt(Math.pow(pt2.x-x, 2) + Math.pow(pt2.y-y, 2));
	}

	public int cross(pt pt2)
	{
		return x*pt2.y - y*pt2.x;
	}

	public boolean isRightTurn(pt mid, pt next)
	{
		pt v1 = vect(mid);
		pt v2 = mid.vect(next);
		return v1.cross(v2) >= 0;
	}

	// used to get rid of collinear points that are already accounted for by
	// points on either end
	public boolean isStraight(pt mid, pt next)
	{
		pt v1 = vect(mid);
		pt v2 = mid.vect(next);
		return v1.cross(v2) == 0;
	}

	public boolean isZero()
	{
		return x == 0 && y == 0;
	}

	public int compareTo(pt other)
	{
		pt v1 = ref.vect(this);
		pt v2 = ref.vect(other);

		if (v1.isZero()) return -1;
		if (v2.isZero()) return 1;

		if (v1.cross(v2) != 0)
			return -v1.cross(v2);

		if (ref.dist(v1) < ref.dist(v2))
			return -1;
		return 1;
	}
}
