module PhyloTreeLikelihood
(
    likelihood,
    treeLikelihood,
    kimura,
    Nucleotide (..),
    PhyloTree (..),
    sequenceLength
) where

import Debug.Trace

debug = flip trace

data Nucleotide = A | G | C | T deriving(Eq, Show)

isPurine :: Nucleotide -> Bool
isPurine A = True
isPurine G = True
isPurine _ = False

isPyrimidine :: Nucleotide -> Bool
isPyrimidine C = True
isPyrimidine T = True
isPyrimidine _ = False

agct = [A,G,C,T]

-- an evolution model takes the two nucleotides at the ends of a branch and the
-- branch length. it returns the prob. of seeing the first nucleotide at one
-- edge of the branch and the second at the other end, given the branch length
type EvModel = (Nucleotide -> Nucleotide -> Double -> Double)

data PhyloTree = Leaf [Nucleotide]
     | Node {
         left        :: (PhyloTree, Double),
         right       :: (PhyloTree, Double),
         likelihoods :: [[Double]]
       }
       deriving (Show)

kimura :: Double -> EvModel
kimura r s x t
    | isTransition   = 1/4-exp(-(2*r+1)/(r+1)*t)/2+exp(-2/(r+1)*t)/4
    | isTransversion = (1/2-exp(-2/(r+1)*t)/2)/2
    | otherwise      = 0.25+exp(-(2*t)/(r+1))/4+exp(((-2*r-1)*t)/(r+1))/2
    where isTransition = ((isPurine s && isPurine x) || (isPyrimidine s && isPyrimidine x)) && s /= x
          isTransversion = (isPurine s && isPyrimidine x) || (isPyrimidine s && isPurine x)

-- the likelihoods at site i
likelihood :: PhyloTree -> EvModel -> [[Double]]
likelihood (Leaf ns) _ = map (\n -> map (\n' -> if n' == n then 1 else 0) agct) ns
likelihood t@(Node (l,ll) (r,rl) [ls]) p
    | (not . null) ls = [ls]
    | otherwise   = map (\i -> map (\n -> l_k n i) agct) [0..(seqLen-1)]
    where seqLen     = sequenceLength t
          l_k s i    = sum (map (left i s) agct)*sum (map (right i s) agct)  `debug` ("site " ++ (show i))
          left  i s' x = p s' x ll*likelihood l p !!i!!nuc x
          right i s' x = p s' x rl*likelihood r p !!i!!nuc x
          nuc n      = case n of { A -> 0; G -> 1; C -> 2; T -> 3 }

treeLikelihood :: PhyloTree -> EvModel -> Double
treeLikelihood t m = sum $ map (\l -> logBase 10 $ sum l/4) $ likelihood t m

sequenceLength :: PhyloTree -> Int
sequenceLength (Leaf ns)        = length ns
sequenceLength (Node (l,_) _ _) = sequenceLength l
