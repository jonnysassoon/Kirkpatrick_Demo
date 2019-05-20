// p1 -> q1
float crossProduct(Point p1, Point q1, Point p2, Point q2) {
  return (q1.x - p1.x) * (q2.y - p2.y) - (q1.y - p1.y) * (q2.x - p2.x);
}

float dotProduct(Point p1, Point q1, Point p2, Point q2) {
  return (q1.x - p1.x) *(q2.x - p2.x) + (q1.y - p1.y) * (q2.y - p2.y);
}

float distance(Point p1, Point p2) {
  return sqrt( pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2));
}

float getAngle(Point p1, Point vertex, Point p2) {
  float cp = crossProduct(vertex, p1, vertex, p2);
  float yVal = cp / (distance(p1, vertex) * distance(p2, vertex));
  float dp = dotProduct(vertex, p1, vertex, p2);
  float xVal = dp / (distance(p1, vertex) * distance(p2, vertex));
  if (xVal >= 0 && yVal >=0) { // 1st quadrant
    return asin(yVal);
  } else if (xVal >=0 && yVal < 0) { // 4th quadrant
    return asin(yVal);
  } else if (xVal < 0 && yVal >= 0) { // 2nd quadrant
    return acos(xVal);
  }
  return -acos(xVal); // 3rd quadrant
}

Point normalize(Point p) {
  float magnitude = sqrt( p.x * p.x + p.y * p.y);
  p.x /= magnitude;
  p.y /= magnitude;
  return p;
}