using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using System.IO;


namespace Algorithm
{
    class AllGraph
    {
        class A {
            public List<int> x;
            public  string y;
        }
        static void Main(string[] args) 
        {
            //NodeList<int> head = new NodeList<int>(0);
            //head.add(1);
            //head.add(2);
            //while (head != null) {
            //    Console.WriteLine(head.val);
            //    head = head.next;
            //}

            string path = @"C:/Users/Administrator/Desktop/testGraph.txt";
            string name = @"C:/Users/Administrator/Desktop/allPaths2.csv";

            //Graph<int> graph = new Graph<int>(path);
            //DFS s = new DFS(graph, 0);
            //Stack<int> pa = new Stack<int>();
            //for (int v = 4; v != 0; v = s.edgeTo[v]) { pa.Push(v); } pa.Push(0);
            //foreach (int i in pa) { Console.Write(i + "-"); } Console.WriteLine();

            //DiGraph G = new DiGraph(name);
            //AllPaths ap = new AllPaths(G, 189, 190);

            //foreach (ArrayList a in ap.path) {
            //    foreach (int i in a) { 
            //        Console.Write(i+"- "); 
            //    }
            //    Console.WriteLine();
            //}

            //DiGraph diG = new DiGraph(path);
            //DirectedCycle huan = new DirectedCycle(diG);
            //foreach (int i in huan.Cycle()) { Console.Write(i); }
            
            //EdgeWeightDigraph G = new EdgeWeightDigraph(path);
            //NetWorkFlow net = new NetWorkFlow();
            //net.MaxFlow(G, 1, 4, 1);
            //Console.WriteLine(net.maxFlow);
            //for (int a = 0; a < net.AEPath.Count;a++)
            //{
            //    foreach (var i in net.AEPath[a])
            //    {
            //        Console.Write(i+"->");
            //    }
            //    Console.Write("\n");
            //}

            Edge[] e = new Edge[8];
            int[] cost = new int[8] { -1, 4, 3, 2, 2, 5, 1, -3 };
            e[0].o = 1; e[0].d = 2;
            e[1].o = 1; e[1].d = 3;
            e[2].o = 2; e[2].d = 3;
            e[3].o = 2; e[3].d = 5;
            e[4].o = 2; e[4].d = 4;
            e[5].o = 4; e[5].d = 3;
            e[6].o = 4; e[6].d = 2;
            e[7].o = 5; e[7].d = 4;
            for (int i = 0; i < 8; i++) { e[i].cost = cost[i]; }
            Bellman_Ford test = new Bellman_Ford();
            test.V = 5; test.E = 8; 
            test.edge = e; test.resurce = 1;
            if (test.bellman_ford() == 0)
            {
                for (int i = 1; i <= 5; i++)
                {
                    test.print(i);
                }
            }
            Console.WriteLine(1);

            Console.WriteLine("----------------");
            A a = new A(); a.x = new List<int>();
            a.x.Add(0); a.y = "ret"; a.x.Add(1); a.x.Add(2);
            List<int> b = a.x; List<int> c = a.x;
            b.RemoveAt(0); Console.WriteLine("a.x[0]  " + a.x[0]);
            c[0] = 34; Console.WriteLine("a.x[0]  " + a.x[0]);

            //ShortestPath.DijkstraSP sp = new ShortestPath.DijkstraSP(G, 0);
            //IEnumerable<ShortestPath.DirectedEdge> pa = sp.PathTo(4);
            //Console.Write("0 -> 4 " + "({0:F2}):", sp.DistTo(4));
            //foreach (ShortestPath.DirectedEdge e in pa)
            //{

            //    string t = e.toString();
            //    Console.Write(t+" ");
            //}


            ///深度、广度 搜索
            //graph.AddEge(0, 1);
            //Console.WriteLine(graph.outToString(4));
            //BFPaths p = new BFPaths(graph, 0);
            //Console.WriteLine(p.hasPathTo(2));
            //foreach (int Path in p.pathTo(4))
            //{
            //    Console.Write(Path+"—");
            //}

            ///连通变量 CC
            //CC test = new CC(graph);
            //int m = test._count();
            //Console.WriteLine(m + " 部分");
            //List<List<int>> comp = new List<List<int>>();
            //for (int i = 0; i < m; i++)
            //{
            //    comp.Add(new List<int>());
            //}
            //for (int v = 0; v < graph.Vnums(); v++)
            //{
            //    comp[test.ID(v)].Add(v);
            //}
            //for (int j = 0; j < m; j++)
            //{
            //    foreach (int k in comp[j])
            //    {
            //        Console.Write(k + " ");
            //    }
            //    Console.WriteLine();
            //}

            ///检测环
            //Cycle c = new Cycle(graph);
            //Console.WriteLine(c.cycleNums());

        }
    }
      /*public class AdjList<T> 
        {
            private int V;//顶点数量
            private int E;//边的数量

            List<Vertex<T>> items;//顶点集合
            public AdjList() : this(10) { }
            public AdjList(int capacity) //指定邻接表的容量 的构造函数
            {
                items = new List<Vertex<T>>(capacity);
            }
            //查找图中是否包含某项
            //Add a Node
            public void AddV(T item) 
            {
                foreach(Vertex<T> v in items) {
                    if (v.data.Equals(item))
                        throw new Exception(item+"已存在于图中");                    
                }
                items.Add(new Vertex<T>(item));
            }
            //Add an Edge
            public void AddE(T from, T to) 
            {
            
            }


            //链表中的结点
            public class Node 
            {
                public Vertex<T> adjV;//邻接点 域
                public Node next;//引用域
                public Node(Vertex<T> value) 
                {
                    adjV = value;
                }
            }
            //存放在数组中的表头结点
            public class Vertex<T> 
            {
                public T data;
                public Node firstEdge;
                public Boolean visited;//访问结点
                public Vertex(T value) //构造函数
                {
                    data = value;
                }

            }
             
        }*/
      //链表
        public class NodeList<T>:IEnumerable<T>
      {
          //private NodeList<T> first;
          public T val;
          public NodeList<T> next;
          public bool visited=false;//默认设为未访问
          public NodeList(T item) 
          { 
              val = item; 
          }
          public NodeList() { next = null; }
          //public NodeList<T> First{
          //    get { return next; }
          //    set { next = value; }
          //} 
          public void add(T item) 
          {
              NodeList<T> oldNode = next;
              //newNode.val = item;
              next = new NodeList<T>();
              next.val = item;
              next.next = oldNode;
          }

          public IEnumerator<T> GetEnumerator()
          {
            NodeList<T> current = next;
            while (current != null)
            {
                //yield return current;
                 yield return  current.val;
                current = current.next;
            }
          }
          IEnumerator IEnumerable.GetEnumerator() {
              return GetEnumerator();
          }
       }
        

      //邻接表（一个元素为链表的List）
      public class Graph<T> 
      {
          private int V;
          private int E;
          public List<NodeList<int>> adj;
          public  Graph(int V) //适用于手动输入图的点，边，小规模
          {
              this.V = V; this.E = 0;
              adj = new List<NodeList<int>>();//创建顶点数为V的邻接表
              for (int v = 0; v < V; v++)
              {
                  adj.Add(new NodeList<int>(v));//所有链表初始化为空，但这样，每个链表的头结点的val=0;所以初始化为 val = v  
                  //每个链表的第0位是顶点，之后才是其邻接点，这里将顶点与其邻接点存在同一个链表中
              }              
          }

          public Graph(string path)//直接读取数据文件
          {
              StreamReader txt = new StreamReader(path);
              //adj = new List<NodeList<int>>(data[0]);
              this.V = int.Parse(txt.ReadLine());
              int E = int.Parse(txt.ReadLine());
              adj = new List<NodeList<int>>();//创建顶点数为V的邻接表
              for (int j = 0; j < V; j++)
              {
                  adj.Add(new NodeList<int>(j));//所有链表初始化为空，但这样，每个链表的头结点的val=0;所以初始化为 val = v  
                  //每个链表的第0位是顶点，之后才是其邻接点，这里将顶点与其邻接点存在同一个链表中
              }
              string line = txt.ReadLine();
              while (line != null)
              {

                  string[] point = line.Split(' ');
                  int v;
                  int w;

                  for (int i = 0; i < point.Length - 1; i++)
                  {
                      v = int.Parse(point[i]);
                      w = int.Parse(point[++i]);
                      AddEge(v, w);
                  }
                  line = txt.ReadLine();

              }

          }

          public void AddEge(int v, int w) 
          {
              adj[v].add(w); adj[w].add(v); E++;
          }
          public int Vnums() { return V; }
          public int Enums() { return E; }

         
          
          //输出到指定顶点(V-1)为止的图
          public string outToString(int V) 
          {
              string s = V + "顶点" + E + "边\n";
              for (int v = 0; v < V; v++)
              {  //每个链表的第0位是顶点，之后才是其邻接点，这里将顶点与其邻接点存在同一个链表中
                 //所以输出时如下
                  s += this.adj[v].val + ": ";
                  this.adj[v] = this.adj[v].next;
                  while (this.adj[v] != null) 
                  {
                      s += this.adj[v].val + " ";
                      this.adj[v] = this.adj[v].next;
                  }
                  s += "\n";
              }
              return s;
          }
      }

      public class DFPaths 
      {
          private Boolean[] marked;//用来记录顶点是否被访问
          private NodeList<int> S;
          private int[] edgeTo;//从起点到一个顶点的已知路径上的最后一个顶点
          public DFPaths(Graph<int> G, int s) //构造函数
          {
              marked = new Boolean[G.Vnums()];
              S = new NodeList<int>(s);
              edgeTo = new int[G.Vnums()];
              dfs(G, s);
          }

          private void dfs(Graph<int> G, int v) //DFS搜索
          {
              G.adj[v].visited = true; marked[v] = true;
              NodeList<int> w = G.adj[v].next;           
              while (w != null) {
                  if (G.adj[w.val].visited == false) {
                      marked[w.val] = true;
                      edgeTo[w.val] = v;
                      dfs(G, w.val);//递归，隐式用来堆栈，而bfs显示用队列
                  }
                  w = w.next;
              }
          }

          public Boolean hasPathTo(int v) //是否存在 s 到 v 的路径
          {
              return marked[v];
          }

          public List<int> pathTo(int v) //s 到 v 的路径，不存在返回 null
          {
              if (!hasPathTo(v)) return null;
              //Stack<int> path = new Stack<int>();
              //for (int c = v; c != s; c = edgeTo[c]) {
              //    path.Push(c);
              //}
              //path.Push(s);
              List<int> path = new List<int>();
              for (int c = v; c != S.val; c = edgeTo[c])
                  path.Add(c);
              path.Add(S.val);
              path.Reverse();
              return path;

          }

          
      }
      public class BFPaths 
      {
          private Boolean[] marked;
          private NodeList<int> S;
          private int[] edgeTo;
          public BFPaths(Graph<int> G, int s) {
              marked = new Boolean[G.Vnums()];
              S = new NodeList<int>(s);
              edgeTo = new int[G.Vnums()];
              bfs(G, s);
          }
          private void bfs(Graph<int> G, int s) {
              Queue<int> visitedNotchecked = new Queue<int>();
              G.adj[s].visited = true; marked[s] = true;
              visitedNotchecked.Enqueue(s);
              while (visitedNotchecked.Count != 0) {
                  int v = visitedNotchecked.Dequeue();//取队列中下一个顶点；删去
                  //
                  NodeList<int> w = G.adj[v].next;//因为w每次都是取G中的链表，
                  while (w != null) {  //所以下面判断时应用G中对应的点来判断，而非 w.visited(w是新创建的链表对象，对它改变无法影响G中元素)
                      if (G.adj[w.val].visited == false) {
                          G.adj[w.val].visited = true; marked[w.val] = true;
                          edgeTo[w.val] = v;
                          visitedNotchecked.Enqueue(w.val);
                      }
                      w = w.next;
                  }//
              }
          }
          public bool hasPathTo(int v) {
              return marked[v];
          }
          public Stack<int> pathTo(int v) {
              if (!hasPathTo(v)) return null;
              Stack<int> Path = new Stack<int>();
              for (int i = v; i != S.val; i = edgeTo[i]) {
                  Path.Push(i);
              }
              Path.Push(S.val);
              return Path;
          }

      }
      public class CC 
      {
          private bool[] marked;
          private int[] id;//连通分量标识符
          private int count;//连通分量数

          public CC(Graph<int> G) 
          {
              marked = new bool[G.Vnums()];
              id = new int[G.Vnums()];
              for (int s = 0; s < G.Vnums(); s++) 
              {
                  if (!marked[s]) {
                      dfs(G, s);
                      count++;
                  }
              }
          }
          private void dfs(Graph<int> G, int v) {
              marked[v] = true;
              id[v] = count;

              NodeList<int> w = G.adj[v].next;
              while (w != null) {
                  if (!marked[w.val]) {
                      dfs(G, w.val);
                  }
                  w = w.next;
              }
          }

          public bool connected(int v, int w) { return id[v] == id[w];}

          public int ID(int v) { return id[v];}

          public int _count() { return count; }
      }
        //未完
      public class Cyclenums 
      {
          private bool[] marked;
          private int C_count;
          public bool hasCycle;
          public Cyclenums(Graph<int> G) {
              marked = new bool[G.Vnums()];
              for (int s = 0; s < G.Vnums(); s++)
              {
                  if (!marked[s]) 
                      dfs(G, s, s);                                        
                  
              }
          }
          private void dfs(Graph<int> G, int v, int u) {
              marked[v] = true;
              NodeList<int> w = G.adj[v];
              while (w != null) 
              {
                  if (!marked[w.val])
                      dfs(G, w.val, v);
                  //else if (w.val != u)
                  //    hasCycle = true;
                  w = w.next;
              }
          }
          public int cycleNums() { return C_count; }
      }

   ///有向图
      public class DiGraph 
      {
          private int V;
          private int E;
          public List<NodeList<int>> adj;

          public DiGraph(int V) 
          {
              this.V = V; this.E = 0;
              adj = new List<NodeList<int>>();//再次注意 List<T> 的初始化
              for (int v = 0; v < V; v++) {
                  adj.Add(new NodeList<int>(v));
              }
          }
          public DiGraph(string path) 
          {
              StreamReader txt = new StreamReader(path);
              this.V = int.Parse(txt.ReadLine());
              this.E = int.Parse(txt.ReadLine());
              adj = new List<NodeList<int>>();
              for (int v = 0; v < V; v++)
              {
                  adj.Add(new NodeList<int>(v));
              }

              string line = txt.ReadLine();
              while (line != null) 
              {
                  string[] sArry = line.Split(new char[2] { ' ', ',' });

                  int v; int w;
                  for (int s = 0; s < sArry.Length - 1; s++) {
                      v = int.Parse(sArry[s]);
                      w = int.Parse(sArry[++s]);
                      addEdge(v, w);
                  }

                  line = txt.ReadLine();
              }

          }

          public int Vnums() { return V; } public int Enums() { return E; }

          public void addEdge(int v, int w) 
          {
              adj[v].add(w);
              E++;
          }
          public NodeList<int> adjOut(int v)//出边集合（链表）
          {
              return adj[v];
          }

          public DiGraph reverse() //反转，为了得到入边集合
          {
              DiGraph R = new DiGraph(V);
              for (int v = 0; v < V; v++) {
                  NodeList<int> w = adj[v];
                  while (w != null) {
                      R.addEdge(w.val, v);
                      w = w.next;
                  }                  
              }
              return R;
          }
      }

      public class DirectedDFS 
      {
          private bool[] marked;

          public DirectedDFS(DiGraph G, int s) 
          {
              marked = new bool[G.Vnums()];
              dfs(G, s);
          }
          public DirectedDFS(DiGraph G, NodeList<int> sources) //链表sources 中的点 可达的所有点
          {
              marked = new bool[G.Vnums()];
              while (sources != null) 
              {
                  if (!marked[sources.val]) {
                      dfs(G, sources.val);
                  }
                  sources = sources.next;
              }
          }

          private void dfs(DiGraph G, int v) 
          {
              marked[v] = true;
              NodeList<int> w = G.adj[v];
              while (w != null) 
              {
                  if (!marked[w.val]) {
                      marked[w.val] = true;
                      dfs(G, w.val);
                  }
              }
          }
          public bool reachable(int v) { return marked[v]; }
      
      }
      public class DFDirectedPaths { }
      public class BFDirectedPaths { }

        //有向图 中的环
      public class DirectedCycle 
      {
          private bool[] marked;
          private int[] edgeTo;
          private Stack<int> cycle;//有向环中的所有点
          private bool[] onStack;//递归调用的栈 上的所有点

          public DirectedCycle(DiGraph G) 
          {
              marked = new bool[G.Vnums()];
              edgeTo = new int[G.Vnums()];
              onStack = new bool[G.Vnums()];
              for (int v = 0; v < G.Vnums(); v++) {
                  if (!marked[v]) {
                      dfs(G, v);
                  }
              }
          }
          private void dfs(DiGraph G, int v) 
          {
              onStack[v] = true; marked[v] = true;
              foreach (int w in G.adj[v]) 
              {
                  if (this.hasCycle()) { return; }
                  else if (!marked[w]) {
                      edgeTo[w] = v; 
                      dfs(G, w);
                  }
                  else if (onStack[w]) {
                      cycle = new Stack<int>();
                      for (int x = v; x != w; x = edgeTo[x]) {
                          cycle.Push(x);
                      }
                      cycle.Push(w);
                      cycle.Push(v);
                  }
                  onStack[v] = false;
              }
          }

          public bool hasCycle() { return cycle != null; }
          public Stack<int> Cycle() { return cycle; }
      }
        
      //有向图 基于DFS 的顶点排序
      public class DFOrder 
      {
          private bool[] marked;
          private Queue<int> pre;//前序：递归调用之前，将顶点加入队列
          private Queue<int> post;//后序：........后，..........
          private Stack<int> reversePost;//逆后序：.后，将顶点压入栈
          public DFOrder(DiGraph G) 
          {
              pre         = new Queue<int>();
              post        = new Queue<int>();
              reversePost = new Stack<int>();
              marked = new bool[G.Vnums()];
              for (int v = 0; v < G.Vnums(); v++) 
              {
                  if (!marked[v]) dfs(G, v);
              }
          }
          //重载
          public DFOrder(EdgeWeightDigraph G) {
              pre = new Queue<int>();
              post = new Queue<int>();
              reversePost = new Stack<int>();
              marked = new bool[G.Vnums()];
              for (int v = 0; v < G.Vnums(); v++)
              {
                  if (!marked[v]) dfs(G, v);
              }
          }

          private void dfs(DiGraph G, int v) 
          {
              pre.Enqueue(v);
              marked[v] = true;
              foreach (int w in G.adj[v]) 
              {
                  if (!marked[w]) { dfs(G, w); }
                  post.Enqueue(v);
                  reversePost.Push(v);
              }
          }
          private void dfs(EdgeWeightDigraph G, int v)//重载
          {
              pre.Enqueue(v);
              marked[v] = true;
              foreach (DirectedEdge w in G.Adj(v))
              {
                  if (!marked[w.To()]) { dfs(G, w.To()); }
                  post.Enqueue(v);
                  reversePost.Push(v);
              }
          }
          public IEnumerable<int> Pre() { return pre; }
          public IEnumerable<int> Post() { return post; }
          public IEnumerable<int> ReversePost() { return reversePost; }
      }

      //只有有向图无环时，才可进行拓扑排序
      public class Topological 
      {
          private IEnumerable<int> order;//顶点的拓扑顺序
          public Topological(DiGraph G) 
          {
              DirectedCycle cyclefind = new DirectedCycle(G);
              if (!cyclefind.hasCycle()) {
                  //这里是核心，使用了DFOrder 类 中的 逆后序方法
                  //一幅 有向无环图 的拓扑顺序 即为所有顶点的 逆后序排列
                  DFOrder dfs = new DFOrder(G);
                  order = dfs.ReversePost();
              }
          }
          
          public IEnumerable<int> Order() { return order; }
          public bool isDAG() { return order != null; }
      }
    /****************************************/
    //for shortestPath//
      //有向加权图寻环
      public class EdgeWeightedCycle 
      {
          private bool[] marked;
          private DirectedEdge[] edgeTo;
          private Stack<DirectedEdge> cycle;//有向环中的 所有 边
          private bool[] onStack;

          public EdgeWeightedCycle(EdgeWeightDigraph G) 
          {
              marked = new bool[G.Vnums()]; edgeTo = new DirectedEdge[G.Vnums()];
              onStack = new bool[G.Vnums()];
              for (int v = 0; v < G.Vnums(); v++) {
                  if (!marked[v]) { dfs(G, v); }
              }
          }

          private void dfs(EdgeWeightDigraph G, int v) {
              onStack[v] = true; marked[v] = true;
              foreach (DirectedEdge e in G.Adj(v)) {
                   if (this.HasCycle())       { return;
                 } else if (!marked[e.To()]) { 
                      edgeTo[e.To()] = e; 
                      dfs(G, e.To());                   
                 } else if (onStack[e.To()]) { 
                      cycle = new Stack<DirectedEdge>();
                      DirectedEdge x = e;
                      for (; x.From() != e.To(); x = edgeTo[x.From()]) 
                      { cycle.Push(x); }
                      cycle.Push(x);
                      return;//??
                 }
              }
              onStack[v] = false;
          }

          public bool HasCycle() { return cycle != null; }
          public IEnumerable<DirectedEdge> Cycle() { return cycle; }
      }
      //有向加权图 拓扑排序（实际也是无环的）
      public class EWDGraphTopological 
      {
          private IEnumerable<int> order;

          public EWDGraphTopological(EdgeWeightDigraph G) {
              EdgeWeightedCycle cycleFind = new EdgeWeightedCycle(G);
              if (!cycleFind.HasCycle()) {
                  DFOrder dfs = new DFOrder(G);
                  order       = dfs.ReversePost();
              }
          }

          public IEnumerable<int> Order() { return order; }
          public bool IsDAG()             { return order != null; }
      }
    /****************************************/
    
    
    //非递归 DFS
        public class DFS {
          private bool[] marked;
          public int[] edgeTo;
          private Stack<int> stack;

          public DFS(Graph<int> G, int s) {
              marked = new bool[G.Vnums()]; edgeTo = new int[G.Vnums()];
              stack = new Stack<int>();
              dfs(G, s);
          }
          private void dfs(Graph<int> G, int v) {
              stack.Push(v); marked[v] = true;
              foreach (int w in G.adj[v]) {
                  if (!marked[w]) {
                      stack.Push(w);
                      edgeTo[w] = v; marked[w] = true;                     
                  }                                                           
              }
              while (stack.Count != 0) {
                  int top = stack.Pop();
                  foreach (int w in G.adj[top]) {
                      if (!marked[w]) {
                          stack.Push(w); edgeTo[w] = top;
                          marked[w] = true;
                      }
                  }
              }

          }

      }

///ALL PATHS
        public class AllPaths 
      {
          private int S;
          private bool[] marked;//是否访问;是否入栈
         
          private Stack<int> stack;
          public List<ArrayList> path;

          public AllPaths(DiGraph G, int s,int d) 
          {
              this.S = s;
              marked = new bool[G.Vnums()]; //instack = new bool[G.Vnums()];
              stack = new Stack<int>(); path = new List<ArrayList>();
              dfs(G, s,d);
          }
          private void dfs(DiGraph G, int v,int d) 
          {
              if (!marked[v]) {
                  NodeList<int> w = G.adj[v];                  
                  stack.Push(w.val);//possible error？？
                  marked[w.val] = true; 
              }              
              foreach (int node in G.adj[v]) 
              {
                  if (!marked[node])
                  {
                      stack.Push(node);
                      marked[node] = true; 
                      if (node == d)
                      {                          
                          path.Add(new ArrayList());
                          foreach (int i in stack)
                          {                             
                            path[path.Count - 1].Add(i);
                          }
                          marked[stack.Pop()] = false;
                          continue;   
                      }                                            
                      dfs(G, node, d);                     
                  }                   
              }
              marked[stack.Pop()] = false;
          }

      }

        struct Edge
        {
            public int o;
            public int d;
            public int cost;
        };

        class Bellman_Ford
        {
            public int E;
            public int V;
            public int resurce;
            int[] Pre;
            int[] dist;

            public Edge[] edge;
           public  int bellman_ford()
            {
                //edge = new Edge[E];
                dist = new int[V+ 1];
                Pre = new int[V + 1];
                for (int i = 0; i <= V; i++)
                {
                    dist[i] = 100;
                    Pre[i] = -1;
                }
                dist[resurce] = 0;
                int flag = 0;
                for (int i = 0; i < V - 1; i++)
                {
                    flag = 0;
                    for (int j = 0; j < E; j++)
                    {
                        if (dist[edge[j].d] > dist[edge[j].o] + edge[j].cost)
                        {
                            dist[edge[j].d] = dist[edge[j].o] + edge[j].cost;
                            Pre[edge[j].d] = edge[j].o;
                            flag = 1;
                        }
                    }
                    if (flag == 0) { break; }
                }
                return flag;
            }

           public void print(int t) 
            {
                int v = t;
                int[] Pre = this.Pre;
                //Console.Write(t);
                while (v != resurce) {
                    Console.Write(v + "->");
                    v = Pre[v];
                }
                Console.Write(resurce);
                Console.WriteLine();
            }

        }

                         
                                          
    
   
}
