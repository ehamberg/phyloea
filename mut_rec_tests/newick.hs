-- Newick parser is taken from https://github.com/rrnewton/PhyBin

import Text.ParserCombinators.Parsec
import Data.List (intercalate)

type BranchLen = Double
type Label = String

-- Even though the Newick format allows it, ignoring interior node
-- labels. (They are not commonly used.)
data NewickTree = NTLeaf BranchLen Label
                | NTInterior BranchLen [NewickTree]
                  deriving (Show, Eq, Ord)

newickText :: NewickTree -> String
newickText (NTInterior bl ts) = "(" ++ interior ts ++ ");"
  where nt (NTLeaf bl lbl)    = lbl ++ ":" ++ show bl
        nt (NTInterior bl ts) = "(" ++ interior ts ++ "):" ++ show bl
        interior              = intercalate "," . map nt

tag l s =
  case s of
    NTLeaf _ n      -> NTLeaf l n
    NTInterior _ ls -> NTInterior l ls

-- This parser ASSUMES that whitespace has been prefiltered from the input.
newickParser :: Parser NewickTree
newickParser = do x <- subtree
                  l <- len
                  char ';'
                  return $ tag l x

-- a [sub] tree is an internal node or a leaf
subtree :: Parser NewickTree
subtree = internal <|> leaf

-- a leaf without
leaf :: Parser NewickTree
leaf = do n <- name
          return $ NTLeaf 0.0 n

internal :: Parser NewickTree
internal = do char '('
              bs <- branchset
              char ')'
              nm <- name -- IGNORED
              return $ NTInterior 0.0 bs

branchset :: Parser [NewickTree]
branchset =
    do b <- branch <?> "at least one branch"
       rest <- option [] $ try $ char ',' >> branchset
       return (b:rest)

branch :: Parser NewickTree
branch = do s<-subtree
            l<-len;
            return $ tag l s

len :: Parser Double
len = option 0.0 $ do char ':'
                      number

number :: Parser Double
number =
  do sign <- option "" $ string "-"
     fst <- many1 digit
     snd <- option "0" $ try $ do char '.'; many1 digit
     return (read (sign ++ fst++"."++snd) :: Double)

name :: Parser String
name = option "" $ many1 (letter <|> digit <|> oneOf "_.-")

main = let (Right tree) =  parse newickParser "" "((A:0.2,B:0.3):0.4,C);"
           in putStrLn $ newickText tree
