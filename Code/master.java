// Adam Doussan AD844156 04/13/2017

public class CombPerm
{
	// ans needs to be intialized before calling, take.length needs to equal
	// words.length
	public static <T> void genCombs(T[] words, ArrayList<ArrayList<T>> ans,
									boolean [] take, int depth)
	{
		if(depth == words.length)
		{
			ArrayList<T> temp = new ArrayList<>();

			for(int i = 0; i < take.length; i++)
				if(take[i])
					temp.add(words[i]);

			ans.add(temp);
			return;
		}

		take[depth] = true;
		genCombs(words,ans,take,depth+1);

		take[depth] = false;
		genCombs(words,ans,take,depth+1);
	}

	// ans & now needs to be intialized before calling, take.length needs to equal
	// words.length
	public static <T> void genPerms(T[] words, ArrayList<ArrayList<T>> ans,
									boolean [] take, int depth, ArrayList<T> now)
	{
		if(depth == words.length)
		{
			ans.add(new ArrayList<>(now));
			return;
		}

		for(int i = 0; i < words.length; i++)
		{
			if(!take[i])
			{
				take[i] = true;
				now.add(words[i]);
				genPerms(words, ans, take, depth+1, now);
				now.remove(now.size()-1);
				take[i] = false;
			}
		}
	}
}

// Adam Doussan AD844156 04/11/2017

public class Graphs
{
	public static void dfs(ArrayList[] graph, boolean [] visited, int v)
	{
		visited[v] = true;
		for(Integer next : ((ArrayList<Integer>[])graph)[v])
			if(!visited[next])
				dfs(graph, visited, next);
	}

	public static int[] bfs(ArrayList[] graph, int v)
	{
		int n = graph.length;
		int [] distance = new int [n];

		Arrays.fill(distance, -1);
		distance[n] = 0;

		Queue<Integer> q = new ArrayDeque<Integer>();
		q.add(v);

		while(!q.isEmpty())
		{
			int cur = q.remove();
			for(Integer next : ((ArrayList<Integer>[])graph)[cur])
			{
				if(distance[next] == -1)
				{
					distance[next] = distance[cur] + 1;
					q.add(next);
				}
			}
		}

		return distance;
	}
}

// Adam Doussan AD844156 04/13/2017

public class Dijkstras
{
	final public static int oo = (int) 1e9;
	private static int n;

	// pass in a value of true to stop if you want early termination for
	// destination d, false if you want the entire dist array to be filled in
	public static int[] run(ArrayList<Edge>[] g, int s, int d, boolean stop)
	{
		n = g.length;

		boolean [] visited = new boolean[n];
		int [] dist = new int [n];
		Arrays.fill(dist, oo);

		PriorityQueue<Edge> pq = new PriorityQueue<Edge>();
		pq.add(new Edge(s, 0));
		dist[s] = 0;

		while(!pq.isEmpty())
		{
			Edge at = pq.remove();
			if(visited[at.e]) continue;
			visited[at.e] = true;

			if(stop && at.e == d) return dist;

			for(Edge adj : g[at.e])
				if(!visited[adj.e] && adj.w + at.w < dist[adj.e])
					pq.add(new Edge(adj.e, dist[adj.e] = adj.w + at.w));
		}

		return dist;
	}
}

class Edge implements Comparable<Edge>
{
	int e, w;

	public Edge(int e, int w)
	{
		this.e = e; this.w = w;
	}

	public int compareTo(Edge o)
	{
		return this.w - o.w;
	}
}

// Adam Doussan AD844156 04/13/2017

public class BellmanFord
{
	final public static int oo = (int) 1e9;

	// returns null if there is a negative cycle, edges are directional
	public static int[] run(Edge[] eList, int n, int s)
	{
		int [] dist = new int [n];
		Arrays.fill(dist, oo);
		dist[s] = 0;

		for(int i = 0; i < n - 1; i++)
			for(Edge e : eList)
				if(dist[e.v1] + e.w < dist[e.v2])
					dist[e.v2] = dist[e.v1] + e.w;

		for(Edge e : eList)
			if(dist[e.v1] + e.w < dist[e.v2])
				return null;

		return dist;
	}
}

class Edge
{
	public int v1, v2, w;

	public Edge(int v1, int v2, int w)
	{
		this.v1 = v1; this.v2 = v2; this.w = w;
	}
}

// Adam Doussan AD844156 04/13/2017

public class Floyd
{
	public static int [][] path;
	final public static int oo = (int)1e9;

	// gets answer without destroying passed in array
	public static int [][] run(int [][] adj)
	{
		int n = adj.length;
		int [][] m = copy(adj);
		path = new int[adj.length][adj.length];

		for(int i = 0; i < adj.length; i++)
			for(int j = 0; j < adj.length; j++)
				if(adj[i][j] != oo)
					path[i][j] = j;
				else
					path[i][j] = oo;

		for(int k = 0; k < n; k++)
			for(int i = 0; i < n; i++)
				for(int j = 0; j < n; j++)
					if(m[i][k] != oo && m[k][j] != oo)
						if(m[i][k] + m[k][j] < m[i][j])
						{
							m[i][j] = m[i][k] + m[k][j];
							path[i][j] = path[i][k];
						}

		return m;
	}

	public static int [][] copy(int [][] a)
	{
		int [][] res = new int [a.length][a[0].length];
			for(int i = 0; i < a.length; i++)
				for(int j = 0; j < a[0].length; j++)
					res[i][j] = a[i][j];
		return res;
	}

	public static ArrayList<Integer> getPath(int i, int j)
	{
		if(path[i][j] == oo)
			return null;

		ArrayList<Integer> myPath = new ArrayList<>();

		myPath.add(i);

		while(i != j)
		{
			i = path[i][j];
			myPath.add(i);
		}

		return myPath;
	}
}

// Adam Doussan AD844156 05/01/2017

public class myFordFulkerson
{
	public int n, s, t;
	public int [][] g;

	public myFordFulkerson(int n)
	{
		this.n = n+2; this.s = n; this.t = n+1;
		g = new int [n+2][n+2];
	}

	public void addEdge(int v1, int v2, int cap)
	{
		g[v1][v2] = cap;
	}

	public int flow()
	{
		int ans = 0;
		int [] prev = new int [n];

		while(bfs(prev))
		{
			int myNext = prev[t];
			int myFlow = g[myNext][t];

			while(myNext != s)
			{
				myFlow = Math.min(myFlow, g[prev[myNext]][myNext]);
				myNext = prev[myNext];
			}

			ans += myFlow;
			myNext = t;

			while(myNext != s)
			{
				int temp = prev[myNext];
				g[myNext][temp] += myFlow;
				g[temp][myNext] -= myFlow;
				myNext = temp;
			}
		}

		return ans;
	}

	public boolean bfs(int[] prev)
	{
		Queue<Integer> q = new ArrayDeque<>();
		boolean [] visited = new boolean [n];

		q.add(s);
		visited[s] = true;

		while(!q.isEmpty())
		{
			int myNext = q.remove();

			for(int i = 0; i < n; i++)
			{
				if(!visited[i] && g[myNext][i] > 0)
				{
					visited[i] = true;
					q.add(i);
					prev[i] = myNext;
				}
			}
		}

		return visited[t];
	}
}

// Credit Arup Guha

class FordFulkerson
{
	int[][] cap;  boolean[] vis;  int n, s, t, oo = (int)(1E9);

	public FordFulkerson(int size) { n = size + 2;  s = n - 2;  t = n - 1;  cap = new int[n][n]; }
	void add(int v1, int v2, int c) {  cap[v1][v2] = c;  }

	int ff()
	{
		vis = new boolean[n];  int f = 0;
		while (true)
		{
			Arrays.fill(vis, false);
			int res = dfs(s, oo);
			if (res == 0) break;
			f += res;
		}
		return f;
	}

	int dfs(int pos, int min)
	{
		if (pos == t)  return min;
		if (vis[pos])  return 0;
		vis[pos] = true;  int f = 0;

		for (int i = 0; i < n; i++)
		{
			if (cap[pos][i] > 0)  f = dfs(i, Math.min(cap[pos][i], min));
			if (f > 0) { cap[pos][i] -= f;  cap[i][pos] += f;  return f; }
		}
		return 0;
	}
}

// Credit Arup Guha

public class Dinic
{
	ArrayDeque<Integer> q;
	ArrayList<Edge>[] adj;
	int n, s, t, oo = (int)1E9;
	boolean[] blocked;
	int[] dist;

	public Dinic (int N)
	{
		n = N; s = n++; t = n++;
		blocked = new boolean[n];
		dist = new int[n];
		q = new ArrayDeque<Integer>();
		adj = new ArrayList[n];
		for(int i = 0; i < n; ++i)
			adj[i] = new ArrayList<Edge>();
	}

	void add(int v1, int v2, int cap, int flow)
	{
		Edge e = new Edge(v1, v2, cap, flow);
		Edge rev = new Edge(v2, v1, 0, 0);
		adj[v1].add(rev.rev = e);
		adj[v2].add(e.rev = rev);
	}

	boolean bfs()
	{
		q.clear();
		Arrays.fill(dist, -1);
		dist[t] = 0;
		q.add(t);
			
		while(!q.isEmpty()) {
			int node = q.poll();
			if(node == s) 
				return true;
			for(Edge e : adj[node])
			{
				if(e.rev.cap > e.rev.flow && dist[e.v2] == -1)
				{
					dist[e.v2] = dist[node] + 1;
					q.add(e.v2);
				}
			}
		}
		return dist[s] != -1;
	}
		
	int dfs(int pos, int min)
	{
		if(pos == t) 
			return min;
		int flow = 0;
		for(Edge e : adj[pos])
		{
			int cur = 0;
			if(!blocked[e.v2] && dist[e.v2] == dist[pos]-1 && e.cap - e.flow > 0)
			{
				cur = dfs(e.v2, Math.min(min-flow, e.cap - e.flow));
				e.flow += cur;
				e.rev.flow = -e.flow;
				flow += cur;
			}
			if(flow == min)
				return flow;
		}
		blocked[pos] = flow != min;
		return flow;
	}
	
	int flow()
	{
		clear();
		int ret = 0;
		while(bfs())
		{
			Arrays.fill(blocked, false);
			ret += dfs(s, oo);
		}
		return ret;
	}
	
	void clear()
	{
		for(ArrayList<Edge> edges : adj)
			for(Edge e : edges)
				e.flow = 0;
	}
}

class Edge
{
	int v1, v2, cap, flow;
	Edge rev;
	Edge(int V1, int V2, int Cap, int Flow)
	{
		v1 = V1; v2 = V2; cap = Cap; flow = Flow;
	}
}

// Conner 04/29/2017

public class TopologicalSort
{
	public static ArrayList<Integer> run(Graph g)
	{
		int degrees[] = new int[g.numVerts];
		for(int i = 0; i < degrees.length; i++)
		{
			ArrayList<Integer> tempArr = g.matrix[i];

			for(int k = 0; k < tempArr.size(); k++)
				degrees[tempArr.get(k)]++;
		}
	
		Queue<Integer> q = new ArrayDeque<Integer>();

		for(int i = 0; i < g.numVerts; i++)
			if(degrees[i] == 0)
				q.add(i);
		
		int visited = 0;
		ArrayList<Integer> result = new ArrayList<Integer>();

		while(!q.isEmpty())
		{	
			int currentEdge = q.poll();
			result.add(currentEdge);
			
			for(int i = 0; i < g.matrix[currentEdge].size(); i++)
				if(--degrees[g.matrix[currentEdge].get(i)] == 0)
					q.add(g.matrix[currentEdge].get(i));

			visited++;
		}
		
		if(visited != g.numVerts) // cycle
			return null;
		
		return result;
	}

}

class Graph
{
	int numVerts;
	ArrayList<Integer> matrix[];

	public Graph(int numVerts)
	{
		this.numVerts = numVerts;
		matrix = new ArrayList[numVerts];

		for(int i = 0; i < matrix.length; i++)
			matrix[i] = new ArrayList<Integer>();
	}

	public void newEdge(int v1, int v2)
	{
		matrix[v1].add(v2);
	}
}

// Conner  04/29/2017

public class Kruskals
{
	private static int count = 0;
	private static int sort = 0;

	public static Edge[] run(Graph g)
	{
		Arrays.sort(g.edgeArray);
		Edge mst[] = new Edge[g.numVerts];
		for(int i = 0; i < g.numVerts; i++)
		{
			mst[i] = new Edge(0,0,0);
		}
		
		SubTree stArray[] = new SubTree[g.numVerts];
		for(int i = 0; i < g.numVerts; i++)
		{
			stArray[i] = new SubTree();
			stArray[i].parent = i;
			stArray[i].rank = 0;
		}
		
		while(sort < g.numVerts - 1)
		{
			Edge next = g.edgeArray[count++];
			
			int v1Base = g.findSet(stArray, next.v1);
			int v2Base = g.findSet(stArray, next.v2);
			
			if(v1Base != v2Base)
			{
				mst[sort++] = next;
				g.union(stArray, v1Base, v2Base);
			}
		}
		return mst;
	}
}

class Graph
{
	int numVerts;
	int numEdges;
	Edge edgeArray[];
	
	public Graph(int numVerts, int numEdges)
	{
		this.numVerts = numVerts;
		this.numEdges = numEdges;
		edgeArray = new Edge[numEdges];
		for(int i = 0; i < numEdges; i++)
		{
			edgeArray[i] = new Edge(0,0,0);
		}
	}
	public int findSet(SubTree stArray[], int key)
	{
		if(stArray[key].parent != key)
		{
			stArray[key].parent = findSet(stArray, stArray[key].parent);
		}
		
		return stArray[key].parent;
	}
	public void union(SubTree[] stArray, int v1, int v2)
	{
		int v1Base = findSet(stArray, v1);
		int v2Base = findSet(stArray, v2);
		
		if(stArray[v1Base].rank < stArray[v2Base].rank)
			stArray[v1].parent = v2Base;
		
		else if(stArray[v2Base].rank < stArray[v1Base].rank)
			stArray[v2Base].parent = v1Base;
		else
		{
			stArray[v2Base].parent = v1Base;
			stArray[v1Base].rank++;
		}
		
	}
}

class Edge implements Comparable<Edge>
{
		int v1, v2, w;

		public Edge(int v1, int v2, int w)
		{
			this.v1 = v1; this.v2 = v2; this.w = w;
		}

		public int compareTo(Edge e)
		{
			return this.w - e.w;
		}
}

class SubTree
{
	int parent;
	int rank;
}

// Adam Doussan AD844156 04/13/2017

public class MyMath
{
	public static long lcm(long a, long b)
	{
		return a / gcd(a,b) * b;
	}

	public static long gcd(long a, long b)
	{
		return (b == 0) ? a : gcd(b, a%b);
	}

	public static int MCSS(int [] days)
	{
		int maxSum = days[0];
		int runningSum = days[0];

		for(int i = 0; i < days.length; i++)
		{
			if(runningSum < 0)
				runningSum = 0;

			runningSum += days[i];

			if(runningSum > maxSum)
				maxSum = runningSum;
		}

		return maxSum;
	}

	public static boolean[] primeSieve(int n)
	{
		boolean[] isPrime = new boolean[n+1];
		Arrays.fill(isPrime, true);
		isPrime[0] = false;
		isPrime[1] = false;

		for(int i = 2; i*i <= n; i++)
			if(isPrime[i])
				for(int j = i + i; j <= n; j += i)
					isPrime[j] = false;

		return isPrime;
	}

	public static long[] extendedEuclid(long a, long b)
	{
		if(b == 0)
			return new long []{1,0,a};
		else
		{
			long q = a / b;
			long r = a % b;
			long [] rec = extendedEuclid(b,r);
			return new long [] {rec[1], rec[0] - q*rec[1], rec[2]};
		}
	}

	public static BigInteger choose(long x, long y)
	{
		if(y<0 || y>x)
			return BigInteger.ZERO;
		if(y == 0 || y == x)
			return BigInteger.ONE;

		BigInteger ans = BigInteger.ONE;

		for (long i = x - y + 1; i <= x; i++)
		    ans = ans.multiply(BigInteger.valueOf(i));
		for (long j = 1; j <= y; j++)
		    ans = ans.divide(BigInteger.valueOf(j));
		return ans;
	}

	public static long [][] chooseDP(int maxDepth)
	{
		long[][] choose = new long[maxDepth][maxDepth];

		for(int i = 0; i < maxDepth; i++)
		{
			choose[i][0] = choose[i][i] = 1;

			for(int j = 1; j < i; j++) 
				choose[i][j] = choose[i-1][j] + choose[i-1][j-1];
		}

		return choose;
	}

	public static BigInteger sumOfDivisors(ArrayList<pair> primeFac)
	{
		BigInteger sum = new BigInteger("1");

		for(pair a : primeFac)
		{
			BigInteger prime = new BigInteger("" + a.prime);

			BigInteger temp = prime.pow((int)a.exp + 1).subtract(BigInteger.ONE);
			temp = temp.divide(prime.subtract(BigInteger.ONE));

			sum = sum.multiply(temp);
		}

		return sum;
	}

	public static ArrayList<pair> primeFactorize(int n)
	{
		ArrayList<pair> res = new ArrayList<>();
		int div = 2;

		while(div*div <= n)
		{
			int exp = 0;
			while(n % div == 0)
			{
				n /= div;
				exp++;
			}

			if(exp > 0) res.add(new pair(div, exp));
			div++;
		}

		if(n > 1) res.add(new pair(n,1));
		return res;
	}
}

class pair
{
	public int prime, exp;

	public pair(int p, int e)
	{
		prime = p; exp = e;
	}
}

// Credit Arup Guha
// Cleaned up by Adam Doussan AD844156 04/15/2017

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

// Adam Doussan AD844156 04/14/2017

public class KnapSack
{
	// unlim = true means can take unlimited times, otherwise each thing can be
	// taken only once. Add multiple times if multiple limited number of same
	// object. kp[money] contains the answer for money.
	public static long[] solve(int [] cost, int [] value, int money, boolean unlim)
	{
		long [] kp = new long [money+1];
		int n = cost.length;

		if(unlim)
		{
			for(int j = 0; j < n; j++)
			{
				for(int k = cost[j]; k < money+1; k++)
				{
					if(kp[k-cost[j]] + value[j] > kp[k])
						kp[k] = kp[k - cost[j]] + value[j];
				}
			}
		}
		
		else
		{
			for(int j = 0; j < n; j++)
			{
				for(int k = money+1; k >= cost[j] ; k--)
				{
					if(kp[k-cost[j]] + value[j] > kp[k])
						kp[k] = kp[k - cost[j]] + value[j];
				}
			}
		}

		return kp;
	}
}

// Ross 04/30/2017

public class MatrixChain
{
	public static void main(String [] args)
	{
		//Matrix A is 2x4, Matrix B is 4x2, Matrix C is 2x3 and so on.
		int matr[] = new int[] {2, 4, 2, 3, 1, 4};
		int size = matr.length;

		int res = MCM(matr, size);
		System.out.println(res);
	}

	static int MCM(int[] matr, int size)
	{
		int a;
		int memo[][] = new int[size][size];

		for(int i=0; i<size; i++)
		{
			memo[i][i] = 0;
		}

		for(int sub=2; sub <= (size-1); sub++)
		{
			for(int j=1; j <= ((size-1) - sub + 1); j++)
			{
				memo[j][sub+j-1] = Integer.MAX_VALUE;
				for(int k=j; k <= (sub+j-2); k++)
				{
					a = memo[j][k] + memo[k+1][j+sub-1] + (matr[j-1] * matr[k] * matr[sub+j-1]);
					if (a < memo[j][j+sub-1])
					{
						memo[j][j+sub-1] = a;
					}
				}
			}
		}

		return memo[1][size-1];
	}
}

// Ross 04/30/2017

public class LCS
{
	static int LCS(String str1, String str2, int len1, int len2)
	{
		int[][] memo = new int[len1+1][len2+1];

		for(int i=1; i<=len1; i++)
		{
			for(int j=1; j<=len2; j++)
			{
				if((i == 0) || (j == 0))
				{
					memo[i][j] = 0;
				}
				else if(str1.charAt(i-1) == str2.charAt(j-1))
				{
					memo[i][j] = 1 + memo[i-1][j-1];
				}
				else
				{
					int sub = Math.max(memo[i][j-1], memo[i-1][j]);
					memo[i][j] = sub;
				}
			}
		}
		return memo[len1][len2];
	}
}
