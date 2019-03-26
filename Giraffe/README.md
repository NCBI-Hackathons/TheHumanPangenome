
## Specific aim:
A Faster, Better Short-Read Mapper with Hit Chaining.

#### Fast Clustering (Stretch Goal #1)
When creating an alignment from numerous graph pathways, the speed of alignment is slowed exponentially with each node present in "snarl" regions.  These "snarl" regions constitute a bubble on the graph where multiple nodes can be chosen.  Each of these nodes can exponentially increase the number of paths within the "snarl", and thus the alignment time.  We attempt to improve upon this by excluding nodes not associated with a given haplotype.

Our algorithm speeds up the clustering after the search hits are gathered.

#### Haplotype-based Hit Joining (Stretch Goal #2)

## Slides
![Image00](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/00.png)
![Image01](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/01.png)
![Image02](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/02.png)
![Image03](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/03.png)
![Image04](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/04.png)

## Stretch Goals:
* Fast Clustering
* Haplotype-based Hit Joining

## Fallback Goal:
Dump hit coverage by node.
