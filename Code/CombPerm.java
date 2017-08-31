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
