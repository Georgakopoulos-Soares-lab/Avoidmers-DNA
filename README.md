# Avoidmers

What is an avoidmer?

Avoidmers as a concept from DNA originates from a branch of computer science usually referred to as <i> stringology </i>.

Avoidmer is any k-mer that does escape the (eventually) unavoidable Zimin pattern for n=3, i.e.

$$Z_{3} = abacaba.$$

In other words, there is not a homomorphism from ${a, b, c}$ to the nucleic acid alphabet such that $Z_{3}# is homomorphically embedded
within the avoidmer. In other words, there is no factor w of k such that $f(Z_{3}) = w$, for any mapping $f$.

## Extraction

We scan the genome from left to right, attempting a one-base pair extension of the avoidmer sequence. Initially, all small k-mers will be avoidmers. 
We extend until the first base pair that breaks the avoidability property, and we backtrack one base-pair and store the maximum avoidmer, if it was
above the given threshold, i.e. in our case at least 50 base pairs long. Our goal is to find the maximum Zimin avoidmers.

## UCSC Genome Browser Track

Below we provide the Avoidmers track in UCSC Genome-Browser.

In particular, we provide the following tracks:

- A bed file containing all avoidmers of at least 50bp long.
- A track with the density of avoidmers per 1kB size window.
- A track with the density of avoidmers per 500bp size window.
- A track with the density of avoidmers of at least 60bp per 500bp size window.

[uscs-browser](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hub_3671779_hs1&lastVirtModeType=default[…]r1%3A1%2D610000&hgsid=2510082887_wZFzMPKwpFTHuq0psrdqPssw64ay)
