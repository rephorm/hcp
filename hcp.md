This is a response to a [question on Stack Overflow](http://stackoverflow.com/questions/9982458/creating-a-sphere-packing-nearest-neighbor-list-from-integers) about index nearest neighbors in a hexagonal close-packed (HCP) lattice.

*Caveat: Unlike the cubic lattice, there isn't a single obvious way to truncate the hexagonal close packed (HCP) structure. I've chosen a specific way of truncation that makes indexing a bit more straightforward. A different truncation scheme will lead to different indexing, and thus different relationships between nearest neighbor indices.*

# Two dimensions

Before tackling the full 3d lattice, lets look at a single 2d sheet containing 25 'atoms' laid out as follows

![Figure 1][fig1]

Notice that alternating rows are shifted back and forth with respect to eachother. Here is one way we can number them. 

![FIG 2][fig2]

In this example, we have 5 rows with 5 atoms each. To be general, lets call the number of rows H, and the number of "columns" (atoms per row) W.
Then, the index can be calculated in the same fashion as for a square lattice


    index = W * row + col


![FIG 3][fig3]

Likewise, we can get back from index to row and column using truncating division and the modulo operator


    row = index / W
    col = index % W


## Neighbors

Now, lets look at the nearest neighbors. First, ignoring edges, we see that each atom has six nearest neighbors. For atom i, two of the neighbors are easy, namely i-1 and i+1. We can also notice that i+W and i-W are always neighbors of i. The last two are a little trickier, since they depend on the row i is in. Consider atom 7. It is neighbors with atoms 1 and 11, which are i-W-1 and i+W-1. But, atom 12 is neighbors with 8 and 18 (i-W+1 and i+W+1). For atoms in an odd row, we have to subtract one, while for atoms in an even row we need to add one. So, for atoms *in the bulk*, we would have


    def neighbors2d(i, W):
      row = i / W
      return [
              i+1,
              i-1,
              i+W,
              i-W,
              i+W+(-1)**row,
              i-W+(-1)**row,
             ]


After fully accounting for edges, we have


    def neighbors2d(i, W, H):
      row = i / W
      col = i % W

      r = -1 if row % 2 else 1 # r = (-1)**row

      nbors = []

      if col != W-1:
        nbors.append(i+1)
      if col != 0:
        nbors.append(i-1)
      if row != H-1:
        nbors.append(i+W)
      if row != 0:
        nbors.append(i-W)
      if (col != 0 or r > 0) and (col != W-1 or r < 0):
        if row != H-1:
          nbors.append(i+W+r)
        if row != 0:
          nbors.append(i-W+r)

      return nbors


# Three dimensions

When we add a second sheet, we have a choice of how to offset it. I chose the following, which is most consistent with how we offset rows.

![FIG 4][fig4]

By numbering the planes (starting with 0), we again have a simple relation between atom index and (col,row,plane) coordinates.


    index = plane * W * H + row * W + col
    plane = index / (W*H)
    plane_index = index % (W*H)
    row = plane_index / W
    col = plane_index % W


![FIG 5][fig5]

Similar to how we had to account for the alternating shift in rows above, we now have to account for the alternating shift in sheets. But, this can be done in a very similar fashion. I found that it helped to reduce things to `plane_index` to figure out the relationships between neighbors, and how they depended on whether the row and plane were even or odd. The following images demonstrate this.

![Figure 6][fig6b]

Finally, after again accounting for edges (which is a bit more complex this time), we have


    def neighbors(i, W, H, D):
      A = W * H

      plane = i / A
      plane_index = i % A
      row = plane_index / W
      col = plane_index % W

      r = -1 if row % 2 else 1   # (-1)**row
      p = -1 if plane % 2 else 1 # (-1)**plane

      nbors = []

      # first include neighbors in same plane
      if col != W-1: nbors.append(i+1)
      if col != 0:   nbors.append(i-1)
      if row != H-1: nbors.append(i+W)
      if row != 0:   nbors.append(i-W)
      if (col != 0 or r > 0) and (col != W-1 or r < 0):
        if row != H-1: nbors.append(i+W+r)
        if row != 0:   nbors.append(i-W+r)

      # now add neighbors from other planes
      if plane != D-1: nbors.append(i+A)
      if plane != 0:   nbors.append(i-A)

      if (col != 0 or p < 0) and (col != W-1 or p > 0):
        if plane != D-1: nbors.append(i+A-p)
        if plane != 0:   nbors.append(i-A-p)

      if ((col != W - 1 or p > 0 or r < 0) and
          (col != 0 or p < 0 or r > 0) and
          (row != H-1 or p < 0) and
          (row != 0 or p > 0)):
        if plane != D-1:
          nbors.append(i + A + p*W + (r-p)/2) #10
        if plane != 0:
          nbors.append(i - A + p*W + (r-p)/2) #11

      return nbors


To make sure I got the logic correct, I used the following test while writing the above function


    def test_neighbors():
      n = lambda i: set(neighbors(i, 5, 5, 5))

      # test bottom layer
      assert n(0) == set([1,5,6,25,30])
      assert n(2) == set([1,3,7,8,26,27,32])
      assert n(4) == set([3,9,28,29,34])
      assert n(5) == set([0,6,10,30])
      assert n(9) == set([3,4,8,13,14,33,34,38])
      assert n(20) == set([15,16,21,45])
      assert n(21) == set([16,17,20,22,45,46])
      assert n(24) == set([19,23,48,49])

      # test second layer
      assert n(25) == set([0,1,26,30,31,50,51])
      assert n(34) == set([4,9,28,29,33,38,39,54,59])
      assert n(36) == set([7,11,12,31,32,35,37,41,42,57,61,62])
      assert n(49) == set([24,44,48,74])


Note that the test doesn't cover all of the unique types of sites, so, there may still be a corner case somewhere that is wrong. I'll leave checking that as an exercise to the reader.

[fig1]: http://i.imgur.com/tRXUR.png
[fig2]: http://i.imgur.com/gymIl.png
[fig3]: http://i.imgur.com/6phKf.png
[fig4]: http://i.imgur.com/Cos70.png
[fig5]: http://i.imgur.com/zi2oi.png
[fig6]: http://i.imgur.com/r9MLs.png
[fig6b]: http://i.imgur.com/foaLV.png
[fig7]: http://i.imgur.com/tFu47.png
