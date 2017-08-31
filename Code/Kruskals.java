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
