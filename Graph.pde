public class MapTriangle extends Triangle {
  HashSet<MapTriangle> children;
  int depth, lvlRank;
  
  public MapTriangle(Point [] verts, int num, int idn) {
     super(verts, num, idn);
     children = new HashSet<MapTriangle>();
     depth = 0;
  }
  public void addChild(MapTriangle child) {
    if (depth < child.depth+1) depth = child.depth+1;
    children.add(child); 
  }
  public void display(int numTabs) {
    for (int i = 0; i < numTabs; i++) print("\t");
    super.display();
    for (MapTriangle child : children) child.display(numTabs+1);
  }
  public MapTriangle getNextTri(Point p) {
    if (children.size() == 0) return this;
    for (MapTriangle mt : children) {
      if (isInsideTriangle(p, mt.verts[0], mt.verts[1], mt.verts[2])) return mt;
    }
    return null;
  }
  
  public MapTriangle getTri(Point p) { // the point will be inside the vertices of this object
    if (children.size() == 0) return this; // we don't need to keep searching
    //MapTriangle child = null;
    for (MapTriangle mt : children) {
      if (isInsideTriangle(p, mt.verts[0], mt.verts[1], mt.verts[2])) return mt.getTri(p);
    }
    return null;
  }
}