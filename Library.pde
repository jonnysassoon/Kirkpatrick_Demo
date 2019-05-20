public class Point {
  float x, y;
  int id;
  Point (float xPos, float yPos, int label) {
    x = xPos;
    y = yPos;
    id = label;
  }
  Point (float xPos, float yPos) {
    x = xPos;
    y = yPos;
    
  }
  public void display() {
    println("(" + x + ", " + y + ")"); 
  }
  public boolean equals(Point other) {
   return x == other.x && y == other.y; 
  }
}

public class Triangle {
  Point [] verts;
  int faceNum, id;
  Triangle(Point [] v, int fn, int idNum) {
    verts = v;
    faceNum = fn;
    id = idNum;
  }
  public void display() {
    print(faceNum + ": [");
    for (Point p : verts) {
      print(p.id + " "); 
    }
    println("]");
  }
}

// Checks if a given point is in the bounds of a triangle
// Quick code from stack overflow 
// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
boolean isInsideTriangle(Point p, Point left, Point top, Point right) {
  float d1, d2, d3;
  boolean has_neg, has_pos;

  d1 = sign(p, left, top);
  d2 = sign(p, top, right);
  d3 = sign(p, right, left);

  has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
  has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

  return !(has_neg && has_pos);
}

float sign (Point p1, Point p2, Point p3) {
  return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

// Checks if we can add a new point to a polygon
boolean noIntersects(Point newP, ArrayList<Point> currPoly, ArrayList<ArrayList<Point>> otherPolys) {
  if (currPoly.size() == 0) return true;
  Point p1 = currPoly.get(currPoly.size()-1);
  if ( p1.equals(newP) ) return false;
  if (segmentIntersectsWithPolygon(p1, newP, currPoly, false, true)) return false;
  for (ArrayList<Point> poly : otherPolys) {
    if (segmentIntersectsWithPolygon(p1, newP, poly, true, false)) return false;
  }
  return true;
}

// TODO: what if I have a star and I try to "split" from one point directly
// to a neighboring point through the outside of the star?
int insidePolygon(Point p, ArrayList<ArrayList<Point>> polygons) {
 for (int i = 0; i < polygons.size(); i++) {
   if (isInsidePolygon(p, polygons.get(i))) return i; 
 }
 return -1; // not inside a polygon
}

// Checks all line segments for overlaps (likely incredibly inefficient)
// check line p1q1 against all points of a polygon
boolean segmentIntersectsWithPolygon(Point p1, Point q1, ArrayList<Point> polygon, boolean alreadyClosed, boolean samePolygon) {
  if  (!finished && polygon.size() <= 1) return false;
  Point p2, q2;
  for (int i = 0; i < polygon.size()-1; i++) { // all segments besides the closing segment
    p2 = polygon.get(i);
    q2 = polygon.get(i+1);
    if (segmentsIntersect(p1, q1, p2, q2, samePolygon)) return true;
  }
  p2 = polygon.get(polygon.size()-1);
  q2 = polygon.get(0);
  if (alreadyClosed) return segmentsIntersect(p1, q1, p2, q2, samePolygon); // the closing segment
  return false;
}

// algorithm from https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
boolean segmentsIntersect(Point p1, Point q1, Point p2, Point q2, boolean samePolygon) {
  /* no math required point checks */
  // complete overlap -> always illegal
  if ( (p1.equals(p2) && q1.equals(q2)) || (p1.equals(q2) && q1.equals(p2)) ) return true;
  if ( p1.equals(q2) ) return false; // start of one is the end of another, alwas fine
  if ( samePolygon ) { // you're checking against another segment in your own polygon
    if ( completedFirstPolygon ) { // we're on some split so we can't ever connect to ourselves, even the start point
      // if the endpoint is the same as either points, it's an illegal segment
      if ( q1.equals(p2) || q1.equals(q2) ) return true; 
    } else { // we're still checking the outer polygon against itself
      // the assumption is that we cannot click to points on the outer polygon while
      // we are building it. If the condition below is true, it most be a closing segment
      if ( q1.equals(p2) ) return false; 
    }
  }
  if ( !samePolygon ) {
    // checking the segment against some segment in another polygon
    if ( p1.equals(p2) || q1.equals(p2) || q1.equals(q2) ) return false;
  }
  
  // This portion from https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
  int o1 = orientation(p1, q1, p2); 
  int o2 = orientation(p1, q1, q2); 
  int o3 = orientation(p2, q2, p1); 
  int o4 = orientation(p2, q2, q1); 

  // General case 
  if (o1 != o2 && o3 != o4) 
      return true; 

  // Special Cases 
  // p1, q1 and p2 are colinear and p2 lies on segment p1q1 
  if (o1 == 0 && onSegment(p1, p2, q1)) return true; 

  // p1, q1 and q2 are colinear and q2 lies on segment p1q1 
  if (o2 == 0 && onSegment(p1, q2, q1)) return true; 

  // p2, q2 and p1 are colinear and p1 lies on segment p2q2 
  if (o3 == 0 && onSegment(p2, p1, q2)) return true; 

   // p2, q2 and q1 are colinear and q1 lies on segment p2q2 
  if (o4 == 0 && onSegment(p2, q1, q2)) return true; 

  return false; // Doesn't fall in any of the above cases
}

int orientation(Point p, Point q, Point r) { 
  // See https://www.geeksforgeeks.org/orientation-3-ordered-points/ 
  // for details of below formula. 
  float val = (q.y - p.y) * (r.x - q.x) - 
            (q.x - p.x) * (r.y - q.y); 

  if (val == 0) return 0;  // colinear 

  return (val > 0)? 1: 2; // clock or counterclock wise 
}

boolean onSegment(Point p, Point q, Point r) {
  if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) && 
    q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y)) return true; 

  return false; 
} 

boolean orientationIsClockwise(ArrayList<Point> polygon) {
  int minInd = getMinXindex(polygon);
  int prevInd = (minInd == 0) ? polygon.size()-1 : minInd-1;
  int nextInd = (minInd == polygon.size()-1) ? 0 : minInd + 1;
  return isLeftTurn(polygon.get(prevInd), polygon.get(minInd), polygon.get(nextInd)); 
}

int getMinXindex(ArrayList<Point> polygon) {
  int minInd = -1;
  float minXcoord = 999999; // this is way outside the bounds
  for (int i = 0; i < polygon.size(); i++) {
    if (polygon.get(i).x < minXcoord) {
      minInd = i;
      minXcoord = polygon.get(i).x;
    }
  }
  return minInd;
}

boolean isLeftTurn(Point p1, Point p2, Point p3) {
 return crossProduct(p1, p2, p2, p3) > 0; 
}

int getMaxXindex(ArrayList<Point> polygon) {
  int maxInd = -1;
  float maxXcoord = -999999; // this is way outside the bounds
  for (int i = 0; i < polygon.size(); i++) {
    if (polygon.get(i).x > maxXcoord) {
      maxInd = i;
      maxXcoord = polygon.get(i).x;
    }
  }
  return maxInd;
}

int getMinYindex(ArrayList<Point> polygon) {
  int minInd = -1;
  float minYcoord = 999999; // this is way outside the bounds
  for (int i = 0; i < polygon.size(); i++) {
    if (polygon.get(i).y < minYcoord) {
      minInd = i;
      minYcoord = polygon.get(i).y;
    }
  }
  return minInd;
}

int getMaxYindex(ArrayList<Point> polygon) {
  int maxInd = -1;
  float maxYcoord = -999999; // this is way outside the bounds
  for (int i = 0; i < polygon.size(); i++) {
    if (polygon.get(i).y > maxYcoord) {
      maxInd = i;
      maxYcoord = polygon.get(i).y;
    }
  }
  return maxInd;
}

// from: https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
boolean isInsidePolygon(Point p, ArrayList<Point> polygon) {
  // There must be at least 3 vertices in polygon[] 
    if (polygon.size() < 3)  return false; 
  
    // Create a point for line segment from p to infinite 
    Point extreme = new Point((int) 1e6, p.y); 
  
    // Count intersections of the above line with sides of polygon 
    int count = 0, i = 0; 
    do
    { 
        int next = (i+1)%polygon.size(); 
  
        // Check if the line segment from 'p' to 'extreme' intersects 
        // with the line segment from 'polygon[i]' to 'polygon[next]' 
        if (segmentsIntersect(polygon.get(i), polygon.get(next), p, extreme, false)) 
        { 
            // If the point 'p' is colinear with line segment 'i-next', 
            // then check if it lies on segment. If it lies, return true, 
            // otherwise false 
            if (orientation(polygon.get(i), p, polygon.get(next)) == 0) 
               return onSegment(polygon.get(i), p, polygon.get(next)); 
  
            count++; 
        } 
        i = next; 
    } while (i != 0); 
  
    // Return true if count is odd, false otherwise 
    return count%2 == 1;
}

// Returns point on the polygon if it was appropriately
Point snapToPoint(Point p, ArrayList<Point> polygon) {
  for (Point aPoint : polygon) {
    if (distance(p, aPoint) < 10) return aPoint; // if you clicked within the radius of a point
  }
  return null;
}

boolean canCutAtVertex(ArrayList<Point> polygon, int earInd, boolean clockwise) {
  Point prev = polygon.get((earInd == 0) ? polygon.size()-1 : earInd - 1);
  Point curr = polygon.get(earInd);
  Point next = polygon.get( (earInd + 1) % polygon.size() );
  if (isReflex(prev, curr, next, clockwise) ) return false; // can't be a reflex vertex
  for (int i = 0; i < polygon.size() - 3; i++) { // check if any other points are inside this triangle (skip the three being used for the ear)
    Point test = polygon.get( (i + 2 + earInd) % polygon.size() );
    if ( isInsideTriangle(test, prev, curr, next) ) return false; // if there exist one, this can't be an ear
  }
  return true; // if there doesn't exist one, it is an ear
}

boolean isReflex(Point p1, Point p2, Point p3, boolean clockwise) {
  if ( isLeftTurn(p1, p2, p3) ) return !clockwise;
  return clockwise;
}

boolean angleLT(Point center, Point pivot, Point p1, Point p2) {
  float ang1 = getAngle(pivot, center, p1);
  float ang2 = getAngle(pivot, center, p2);
  if (ang1 > 0 && ang2 > 0 || ang1 < 0 && ang2 < 0) {
    return ang1 < ang2;
  } else {
    return ang1 > 0;
  }
}

void swap(ArrayList<Point> arr, int i, int j) {
  Point tmp = arr.get(i);
  arr.set(i, arr.get(j));
  arr.set(j, tmp);
}

void radialSort(ArrayList<Point> neighbors, Point center) { // something wrong with this
  // because neighbors.size() = O(1), do insertion sort
  Point pivot = neighbors.get(0);
  for (int i = 1; i < neighbors.size(); i++) {
    Point curr = neighbors.get(i);
    for (int j = i-1; j > 0 && angleLT(center, pivot, curr, neighbors.get(j)); j--) {
      swap(neighbors, j+1, j); 
    }
  }
}

boolean trianglesOverlap(MapTriangle t1, MapTriangle t2) {
  // either one point is inside the other triangle
  if (isInsideTriangle(t1.verts[0], t2.verts[0], t2.verts[1], t2.verts[2]) ||
      isInsideTriangle(t1.verts[1], t2.verts[0], t2.verts[1], t2.verts[2]) ||
      isInsideTriangle(t1.verts[2], t2.verts[0], t2.verts[1], t2.verts[2])) return true;
  ArrayList<Point> tri2 = new ArrayList<Point>();
  tri2.add(t2.verts[0]);
  tri2.add(t2.verts[1]);
  tri2.add(t2.verts[2]);
  // or some segment overlaps with it
  if ( segmentIntersectsWithPolygon(t1.verts[0], t1.verts[1], tri2, true, false) ||
      segmentIntersectsWithPolygon(t1.verts[1], t1.verts[2], tri2, true, false) ||
      segmentIntersectsWithPolygon(t1.verts[2], t1.verts[0], tri2, true, false)) return true;
  return false;
}


void dispAllPolygons() {
  for (int i = 0; i < completedPolygons.size(); i++) {
   dispPolygon(completedPolygons.get(i)); 
  }
}

void dispPolygon(ArrayList<Point> polygon) {
  print("[");
  for (int i = 0; i < polygon.size(); i++) {
    print(polygon.get(i).id);
    if (i != polygon.size()-1) print(", ");
  }
  println("]");
}

void dispLookup(HashMap<Point, HashMap<String, Object>> lookupTable) {
  Iterator it = lookupTable.entrySet().iterator();
  while (it.hasNext()) {
      Map.Entry pair = (Map.Entry)it.next();
      Point p = (Point) pair.getKey();
      HashMap<String, Object> adj = (HashMap<String, Object>) pair.getValue();
      HashSet<Point> adjNodes = (HashSet<Point>) adj.get("adjNodes");
      HashSet<Triangle> adjTris = (HashSet<Triangle>) adj.get("adjTris");
      print(p.id + ":{\n\tnodes: {");
      for (Point q : adjNodes) {
        print(q.id + " ");
      }
      print("}\n\tTriangles: {");
      for (Triangle tri :adjTris) {
        print("[");
        for (Point q : tri.verts) {
          print(q.id + " ");
        }
        print("] ");
      } println("}");
      println("}");
  } 
}