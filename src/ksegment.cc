#include <config.h>

#include <ksegment.h>

//  K_SEGMENT :: K_SEGMENT()
//    constructs a null segment.

K_SEGMENT :: K_SEGMENT()
{
  start = 0;
  end   = 0;

  ref_count = 0;
}

//  K_SEGMENT :: K_SEGMENT(const K_POINT2D& x, const K_POINT2D& y)
//    constructs a segment (x, y].

K_SEGMENT :: K_SEGMENT(const K_POINT2D& x, const K_POINT2D& y)
{
  start = new K_POINT2D(x);
  start->ref_count++;
  end   = new K_POINT2D(y);
  end->ref_count++;

  ref_count = 0;
}

//  K_SEGMENT :: K_SEGMENT(K_POINT2D* const x, K_POINT2D* const y)
//    constructs a segment (x, y].

K_SEGMENT :: K_SEGMENT(K_POINT2D* const x, K_POINT2D* const y)
{
  start = x;
  start->ref_count++;
  end   = y;
  end->ref_count++;

  ref_count = 0;
}

//  K_SEGMENT :: K_SEGMENT(const K_SEGMENT& s)
//    the copy constructor

K_SEGMENT :: K_SEGMENT(const K_SEGMENT& s)
{
  start = s.start;
  start->ref_count++;
  end   = s.end;
  end->ref_count++;

  ref_count = 0;
}

//  K_SEGMENT& K_SEGMENT :: operator =(const K_SEGMENT& s)
//    the assignment operator

K_SEGMENT& K_SEGMENT :: operator =(const K_SEGMENT& s)
{
  if (this != &s)
  {
    if (start && !--start->ref_count)
      delete start;

    if (end && !--end->ref_count)
      delete end;

    if (start = s.start)
      start->ref_count++;

    if (end = s.end)
      end->ref_count++;
  }

  return *this;
}

//  K_SEGMENT :: ~K_SEGMENT()
//    the destructor

K_SEGMENT :: ~K_SEGMENT()
{
  if (start && !--start->ref_count)
    delete start;

  if (end && !--end->ref_count)
    delete end;
}

ostream& K_SEGMENT :: output(ostream& o) const
{
  if (start && end)
    o << "[ " << *start << ", " << *end << " )" << flush;
  else  //  if (!start || !end)
    o << " NULL " << flush;

  return o;
}

ostream& operator <<(ostream& o, const K_SEGMENT& s)
{
  return s.output(o);
}

//  K_SEGMENT K_SEGMENT :: reverse() const
//    return s = (*end, *start].

K_SEGMENT K_SEGMENT :: reverse() const
{
  K_SEGMENT s;

  if (s.start = end)
    s.start->ref_count++;

  if (s.end = start)
    s.end->ref_count++;

  return s;
}

//  K_BOXCO2 K_SEGMENT :: outer_box() const
//    returns the outer box for *this.

K_BOXCO2 K_SEGMENT :: outer_box() const
{
  unsigned long i;
  K_BOXCO2      b;

  b = start->bbox().merge(end->bbox());

  for (i = 0; i < 2; i++)
    if (b.low[i] == b.high[i])
      b.low_open[i] = b.high_open[i] = 0;
    else  //  if (b.low[i] != b.high[i])
      b.low_open[i] = b.high_open[i] = 1;

  return b;
}

//  K_BOXCO2* K_SEGMENT :: inner_box() const
//    returns a pointer to the inner box for *this.
//    returns 0 if the inner box for *this is not well-def'ed.

K_BOXCO2* K_SEGMENT :: inner_box() const
{
  unsigned long i;
  int           well_defed;
  K_BOXCO2      s, e, b_proto;
  K_BOXCO2*     b;

  s          = start->bbox();
  e          = end->bbox();
  well_defed = 1;
  i          = 0;

  while (well_defed && i < 2)
    if (s.low[i] > e.high[i])
    {
      b_proto.low[i]       = e.high[i];
      b_proto.high[i]      = s.low[i];
      b_proto.low_open[i]  = 1 - e.high_open[i];
      b_proto.high_open[i] = 1 - s.low_open[i];
      i++;
    }
    else if (s.high[i] < e.low[i])
    {
      b_proto.low[i]       = s.high[i];
      b_proto.high[i]      = e.low[i];
      b_proto.low_open[i]  = 1 - s.high_open[i];
      b_proto.high_open[i] = 1 - e.low_open[i];
      i++;
    }
    else if (s.low[i] == e.high[i] && s.low_open[i] && e.high_open[i])
    {
      b_proto.low[i]      = b_proto.high[i]      = s.low[i];
      b_proto.low_open[i] = b_proto.high_open[i] = 0;
      i++;
    }
    else if (s.high[i] == e.low[i] && s.high_open[i] && e.low_open[i])
    {
      b_proto.low[i]      = b_proto.high[i]      = s.high[i];
      b_proto.low_open[i] = b_proto.high_open[i] = 0;
      i++;
    }
    else
//    {
//      cerr << " ksegment: inner_box: inner_box is not well-def'ed! " << endl << flush;
      well_defed = 0;
//    }

  if (well_defed)
    b = new K_BOXCO2(b_proto);
  else  //  if (!well_defed)
    b = 0;

  return b;
}

//  int K_SEGMENT :: contains(K_POINT2D& x)
//    returns 1 if x lies inside the inner box of *this,
//            0 if x lies outside the outer box of *this.
//    POSSIBLY DOES NOT TERMINATE!
//  CAUTION!!!  The start point is OUT of the segment whereas
//              the end point is IN the segment.

int K_SEGMENT :: contains(K_POINT2D& x)
{
  int c, o;

  if (start->equal(x))  //  x is the start point of *this.
    c = 0;
  else  //  if (!start->equal(x))
  {
    start->separate(x);  //  !start->equal(x) => start->separate(x) terminates.

    if (end->equal(x))  //  x is the end point of *this.
      c = 1;
    else  //  if (!end->equal(x))
          //  x is neither the start point nor the end point of *this.
    {
      end->separate(x);  //  !end->equal(x) => end->separate(x) terminates.

      K_BOXCO2 obox, xbox;

      obox = outer_box();
      xbox = x.bbox();

      if ((o = obox.overlap(xbox)) && !obox.contains(xbox))
      //  x overlaps the outer box of *this but is not entirely inside.
      //  Cut x to see if x is inside or outside the outer box of *this.
      {
        x.cut_s(obox.low[0]);
        x.cut_s(obox.high[0]);
        x.cut_t(obox.low[1]);
        x.cut_t(obox.high[1]);
        xbox = x.bbox();

        o = obox.overlap(xbox);
      }

      //  x is either entirely outside the outer box of *this or
      //              entirely inside.

      if (!o)  //  x is entirely outside the outer box of *this.
        c = 0;
      else  //  if (o)
            //  x is entirely inside the outer box of *this.
        if (obox.low[0] == obox.high[0])
        //  x is on *this AND neither the start point nor the end point.
        {
          assert(x.type == 3 || x.type == 4);
          c = 1;
        }
        else if (obox.low[1] == obox.high[1])
        //  x is on *this AND neither the start point nor the end point.
        {
          assert(x.type == 2 || x.type == 4);
          c = 1;
        }
        else if (start->equal_s(*end))
        //  *this is vertical; c = 1 iff (*start, *end) // (*start, x)
          c = start->equal_s(x);
        else if (start->equal_t(*end))
        //  *this is horizontal; c = 1 iff (*start, *end) // (*start, x)
          c = start->equal_t(x);
        else
        //  *this is neither vertical nor horizontal.
        //  Thus, its inner box is well-defined.
        //  Cut x to see x is inside or outside the inner box of *this.
        {
          K_BOXCO2* ibox;

          if (ibox = inner_box())
          {
            x.cut_s(ibox->low[0]);
            x.cut_s(ibox->high[0]);
            x.cut_t(ibox->low[1]);
            x.cut_t(ibox->high[1]);
            xbox = x.bbox();

            if (ibox->contains(xbox))
            //  x is entirely inside the inner box of *this.
              c = 1;
            else
            //  x is inside the outer box but outside the inner box of *this.
              c = - 1;

            delete ibox;
          }
          else
          //  ibox is not well-defined.
            c = - 1;

          //  Shrink start and end until x is either in or out.

          while (c < 0)
          {
//            cerr << " XXX: ksegment: contains: shrinking *start and *end " << endl << flush;
            start->shrink(shrink_step, shrink_step);
            end->shrink(shrink_step, shrink_step);

            if (ibox = inner_box())
            {
              obox = outer_box();
              x.cut_s(ibox->low[0]);
              x.cut_s(ibox->high[0]);
              x.cut_t(ibox->low[1]);
              x.cut_t(ibox->high[1]);
              x.cut_s(obox.low[0]);
              x.cut_s(obox.high[0]);
              x.cut_t(obox.low[1]);
              x.cut_t(obox.high[1]);
              xbox = x.bbox();

              if (ibox->contains(xbox))
              //  x is inside of the inner box of *this.
                c = 1;
              else if (!obox.contains(xbox))
              //  x is outside of the outer box of *this.
                c = 0;
              else
                c = - 1;

              delete ibox;
            }
            else
              c = - 1;
          }
        }
    }
  }

  return c;
}

