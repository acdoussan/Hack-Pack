// Adam Doussan AD844156 04/11/2017

public class Trees
{
	public int height(Node root)
	{
		if(root == null) return -1;
		int res = 0, i;
	
		for(int i = 0; i < root.numChildren; i++)
		{
			int cur = height(root.kids[i]);
			res = Math.max(res, cur+1);
		}

		return res;
	}

	public int maxValue(Node root)
	{
		if(root == null) return 0;

		int leftSide = maxValue(root.left);
		int rightSide = maxValue(root.right);

		if(root.cash > leftSide + rightSide)
			return root.cash;

		return leftSide + rightSide;
	}
}
