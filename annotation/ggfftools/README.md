# ggfftools
> WIP

### Motivation
A common task while working with region based annotations is to perform set operations like intersection, union, difference, etc.
on region lists.

ggfftools tries to implement a few utilities to perform such set operations on gGFF files described [previously](../annotation#the-ggff-format)

### Design
An interval in gGFF format is represented as a set of node slices. [Roaring bitmaps](http://roaringbitmap.org/about/)
are used to store node ids (uint64's) which provide access to fast set operations while maintaining better compression 
ratios than uncompressed or other compressed bitmaps (EWAH).

For simplicity, [leveldb](https://github.com/google/leveldb) is used to store the node slices for each gGFF record.

### Operations
* Set of node ids (or) node slices ---> [QUERY] ---> Set of matching gGFF records in the index
  corresponding to indexed gGFF records matching the given query node slice,
  each record has feature name + gff3 attributes.
  ```
  This currently iterates over all available subgraphs (node slices) in the leveldb store and 
  then performs the required set operation.
  
  But it should be possible using locality sensitive hashing to get a list of subgraphs which are most similar
  to the query node slice and then only perform set operations on those subgraphs.
  ```
  
* Feature name (e.g. "Simple Repeats") ---> [QUERY] ---> Set of matching gGFF records in the index
