// Adam Doussan AD844156 04/29/2017

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
