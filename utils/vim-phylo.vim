""===========================================================================
"" Some Vim tips for working with phylo-data
""
"" Johan Nylander <jnylander@users.sourceforge.net>
""
"" 04/16/2010 12:16:39 PM CEST
""
""
"" Note -- Some functions calls external programs/scripts.
""         Comments/additions/improvements/bug fixes most welcome.
""
"" Discalimer -- Use at your own risk, Caveat emptor!
""
"" See also: http://www.abc.se/~nylander/vim/Align.Vim.tgz
""           http://www.abc.se/~nylander/vim/Nexus.Vim.tgz
""           http://www.abc.se/~nylander/vim/vim-phylo.flv (demo movie in flv format, no sound.)
""


""===========================================================================
"" Count A, C, G, and T's on a line
"" Source -- Alexandru Tudor Constantinescu, 12/14/2004
"" This function counts the AT and GC on the WHOLE line (and only one line)
"" It reports the AT/GC contents and the calculated Tm based on 4*GC+2*AT
"" PLEASE e-mail me (JN) when you have had time to modify this function to count
"" only the "sequence" part of the line.
function! Count_bases() 
   let l:string_length = strlen(substitute(getline("."), ".*", "&", "g"))
   let l:a = l:string_length - strlen(substitute(getline("."), "\\c[a]", "", "g"))
   let l:c = l:string_length - strlen(substitute(getline("."), "\\c[c]", "", "g"))
   let l:g = l:string_length - strlen(substitute(getline("."), "\\c[g]", "", "g"))
   let l:t = l:string_length - strlen(substitute(getline("."), "\\c[t]", "", "g"))
   let l:gap = l:string_length - strlen(substitute(getline("."), "\\c[-]", "", "g"))
   let l:other = l:string_length - l:a - l:c - l:g - l:t - l:gap
   echo "Length=" . l:string_length " (" "A=" . l:a "C=" . l:c "G=" . l:g "T=" . l:t "-=" . l:gap "other=" . l:other ")"
endfunction
"map ZZ :call Count_bases()<CR>


""===========================================================================
"" Align columns based on white space by visually select a text, and then
"" :Align<CR> or press '\a'
"" Source -- ?
"" Good for align sequence labels with unequal length with the seqs.
function! AlignSection(regex) range
  let extra = 1
  let sep = empty(a:regex) ? '\s\+' : a:regex
  let maxpos = 0
  let section = getline(a:firstline, a:lastline)
  for line in section
    let pos = match(line, ' *'.sep)
    if maxpos < pos
      let maxpos = pos
    endif
  endfor
  call map(section, 'AlignLine(v:val, sep, maxpos, extra)')
  call setline(a:firstline, section)
endfunction

function! AlignLine(line, sep, maxpos, extra)
  let m = matchlist(a:line, '\(.\{-}\) \{-}\('.a:sep.'.*\)')
  if empty(m)
    return a:line
  endif
  let spaces = repeat(' ', a:maxpos - strlen(m[1]) + a:extra)
  return m[1] . spaces . m[2]
endfunction
command! -nargs=? -range Align <line1>,<line2>call AlignSection('<args>')
vnoremap <silent> <Leader>a :Align<CR>


""===========================================================================
"" Reverse-Complement DNA sequence
"" Source -- Alexandru Tudor Constantinescu,  02/14/2005
""           Tim Chase and William Nater 12/12/2004
"" 04/13/2010 01:09:14 PM CEST: didn't work as expected
"" This function tries to get the reverse-complement of a certain DNA sequence
"" you have to select a block of text in advance
"" the script is crude, in that it assumes you have only ATCG
"" capitalization does not get screwed up
"" at this moment 12/07/2004, the whole LINE (i.e. not only part of line) gets
"" changed, irrespective of what you select.
"" Tim Chase and William Nater 12/12/2004
"" since ignorecase gives problems (i.e. capitalization is lost)
"" The replacement should be done by selecting a block of text (beware that the
"" WHOLE line will get changed!!) and then issuing the commmand:
"" :RC
"fun! Rev(result)
"   let l:i = strlen(a:result) - 1
"   let l:result = ''
"   while (l:i > -1)
"     let l:result = l:result.a:result[l:i]
"     let l:i = l:i - 1
"   endwhile
"   return l:result
"endfun
"
"function! RC_Tim(l1, l2)
"   let l:str = getline(a:l1)
"   let l:len = strlen(l:str)
"   let l:ignorecs = &l:ic
"   let &l:ic = 0
"   exe a:l1.",".a:l2."j!"
"   exe a:l1."s/.*/\\=Rev(submatch(0))/"
"   exe a:l1."s/\\c[agct]/\\=\"ATGCatgc\"[match(\"TACGtacg\", submatch(0))]/ge"
"   exe a:l1."s/.\\{".&tw."\\}\\zs/\\r/g"
"   let &l:ic = l:ignorecs
"endfunction
"command! -n=* -range RC :call RC_Tim(<line1>,<line2>)


""===========================================================================
"" Nexus-Comment/uncomment
"" Source -- Modified by JN from vim.org/tips
"" Comments/uncomments line by line!
function! NexusComment()
    if getline(".") =~ '['
        let hls=@/
        s/^\(\s*\)\[/\1/
        s/\(\s*\)\]$/\1/
        let @/=hls
    else
        let hls=@/
        s/^/[/
        s/$/]/
        let @/=hls
    endif
endfunction
"map ,[ :call NexusComment()<CR>


""===========================================================================
"" Count the number of '>' in fasta file
"" Source -- JN
function! GetNtax()
    let ntax=0
    g/^/let ntax+=strlen(substitute(getline('.'), '[^>]', '', 'g'))
    return ntax
endfunction

function! GetNtax2()
    exe ': normal G'
    let ntax=0
    g/^/let ntax+=strlen(substitute(getline('.'), '[^>]', '', 'g'))
    echo "Number of sequences: " ntax
endfunction


""===========================================================================
"" Remove gaps in Fasta
"" Source -- JN
function! DegapFasta()
    exe ':g!/^>/s/-//g'
    exe ':g/^\s*$/d'
    exe ':normal gg'
endfunction


""===========================================================================
"" Count the length of longest line in file (not counting white space)
"" Source -- JN from Tim Chase-2, 2007
function! GetMaxLineLength()
    let maxLength=0
    let start=line("1")
    let end=line("$")
    while (start <= end)
        let lineLength=strlen(substitute(getline(start), '\s', '', ''))
        if (lineLength > maxLength)
            let maxLength=lineLength
        endif
        let start = start + 1
    endwhile
    "echo "max line length: " maxLength
    return maxLength
endfunction


""===========================================================================
"" FASTA interleaved to Phyml converter.
"" Source -- JN
function! Fasta2Phyml()
    let ntax=GetNtax()
    exe ':normal gg'
    exe ':s/>/\r>/'
    exe ':%s/\(^>.\+\)$/\1<skojj/'
    exe ':%s/\n//g'
    exe ':%s/>/\r/g'
    exe ':%s/<skojj/<skojj\r/'
    let nchar=GetMaxLineLength()
    exe ':%s/<skojj\n/\t/'
    exe ':normal gg'
    exe ':%Align'
    exe ':normal gg'
    exe "normal i" ntax nchar "\<Esc>"
    "echo "Ntax: " ntax "Nchar: " nchar
endfunction


""===========================================================================
"" Convert interleaved FASTA to non-interleaved FASTA
"" Source -- JN
command! -nargs=0 Fasta2NonInterLeavedFasta :normal gg | :s/>/\r>/<CR> | :%s/\(^>.\+\)$/\1<skojj/<CR> | :%s/\n//g<CR> | :%s/>/\r>/g<CR> | :%s/<skojj/\r/<CR> | :normal ggdd<CR>
command! -nargs=0 UnWrapFasta :normal gg | :Fasta2NonInterLeavedFasta<CR>


""===========================================================================
"" Convert sequential PHYML file format to sequential FASTA file format
"" Source -- JN
command! -nargs=0 Phyml2Fasta :normal ggdd | :%s/^/>/<CR> | :%s/\s\+/\r/<CR>


""===========================================================================
"" Convert FASTA (interleaved) file format to PHYML (sequential) file format
"" Source -- JN
command! -nargs=0 Fasta2Phyml2 :normal gg | :s/>/\r>/<CR> | :%s/\(^>.\+\)$/\1<skojj/<CR> | :%s/\n//g<CR> | :%s/>/\r/g<CR> | :%s/<skojj/\t/<CR> | :normal gg<CR>


""===========================================================================
"" Convert FASTA (non-interleaved) file format to FASTA (interleaved) file format
"" Source -- JN
command! -nargs=0 WrapFasta :normal gg | :g!/^>/s/.\{70}/&\r/g<CR> | :normal gg<CR>


""===========================================================================
"" Convert FASTA file format to sequential PHYML file format
"" Source -- JN
"command! -nargs=0 Fasta2Phyml :normal gg | :%s/^/>/<CR> | :%s/\s\+/\r/<CR>


""===========================================================================
"" FASTA to NEXUS using external program
"" Source -- JN
command! -nargs=0 Fas2Nex :%! fas2nex.stdout %<CR>


""===========================================================================
"" Insert random sequence using external Perl-script
"" Source -- JN
function! RandSeq(len)
    let LENGTH = a:len
    exe ':r!/home/nylander/bin/getrandomsequence.pl' . ' ' . LENGTH
endfunction
command! -nargs=1 Randseq : call RandSeq(<args>)


""===========================================================================
"" Function to read Man page
"" Source -- Edited by JN from http://vim.wikia.com/wiki/
""           Open_a_window_with_the_man_page_for_the_word_under_the_cursor
function! ReadMan(man_word)
    exe ':tabnew'
    exe ':r!man ' . a:man_word . ' | col -b'
    exe ':goto'
    exe ':delete'
    exe ':set filetype=man'
endfunction


""===========================================================================
"" My Phylo Menu (very experimental!) Uses a number of external programs
"" Source -- JN
""===========================================================================
"" Run Alignment programs
menu Phylo.Do\ Alignment.CLUSTALW.Protein :! clustalw -outorder=INPUT -output=GDE -case=UPPER -outfile=ClUsTaL.aln -align -infile=%<CR>: tabe ClUsTaL.aln<CR> :% s/%/>/<CR>: normal gg<CR>
menu Phylo.Do\ Alignment.CLUSTALW.DNA     :! clustalw -outorder=INPUT -output=GDE -case=UPPER -outfile=ClUsTaL.aln -align -infile=%<CR>: tabe ClUsTaL.aln<CR> :% s/#/>/<CR>: normal gg<CR>
"menu Phylo.Do\ Alignment.CLUSTALW.-Sep-       :
menu Phylo.Do\ Alignment.CLUSTALW.Read\ CLUSTALW\ man\ page : call ReadMan('clustalw')<CR>
menu Phylo.Do\ Alignment.MAFFT.mafft       :! mafft  % > MaFfT.mafft.ali<CR> : tabe  MaFfT.mafft.ali<CR><CR>
menu Phylo.Do\ Alignment.MAFFT.linsi       :! linsi  % > MaFfT.linsi.ali<CR> : tabe  MaFfT.linsi.ali<CR><CR>
menu Phylo.Do\ Alignment.MAFFT.ginsi       :! ginsi  % > MaFfT.ginsi.ali<CR> : tabe  MaFfT.ginsi.ali<CR><CR>
menu Phylo.Do\ Alignment.MAFFT.einsi       :! einsi  % > MaFfT.einsi.ali<CR> : tabe  MaFfT.einsi.ali<CR><CR>
menu Phylo.Do\ Alignment.MAFFT.fftnsi      :! fftnsi % > MaFfT.fftnsi.ali<CR>: tabe  MaFfT.fftnsi.ali<CR><CR>
menu Phylo.Do\ Alignment.MAFFT.fftns       :! fftns  % > MaFfT.fftns.ali<CR> : tabe  MaFfT.fftns.ali<CR><CR>
menu Phylo.Do\ Alignment.MAFFT.nwns        :! nwns   % > MaFfT.nwns.ali<CR>  : tabe  MaFfT.nwns.ali<CR><CR>
menu Phylo.Do\ Alignment.MAFFT.nwnsi       :! nwnsi  % > MaFfT.nwnsi.ali<CR> : tabe  MaFfT.nwnsi.ali<CR><CR>
menu Phylo.Do\ Alignment.MAFFT.-Sep-       :
menu Phylo.Do\ Alignment.MAFFT.Read\ MAFFT\ man\ page : call ReadMan('mafft')<CR>

"" Edit alignments
"menu Phylo.Edit.-Sep-  :
menu Phylo.Utilities.Remove\ all\ gaps\ in\ FASTA          : call DegapFasta()<CR>
menu Phylo.Utilities.Align\ taxlabels\ (in\ phyml\ format) : Align<CR>
menu Phylo.Utilities.Count\ sequences\ (in\ FASTA\ format) : echo ''<CR>: call GetNtax2()<CR>
menu Phylo.Utilities.Set\ filetype\ to\ DNA                : set ft=align<CR>
menu Phylo.Utilities.Set\ filetype\ to\ AA                 : set ft=aalign<CR>
menu Phylo.Utilities.Set\ filetype\ to\ Nexus              : set ft=nexus<CR>

"" Format conversions
menu Phylo.Convert.Phyml\ to\ FASTA : Phyml2Fasta<CR>
menu Phylo.Convert.FASTA\ to\ Phyml : call Fasta2Phyml()<CR><CR>
menu Phylo.Convert.FASTA\ to\ Nexus : Fas2Nex<CR> : set ft=nexus<CR>
menu Phylo.Convert.Unwrap\ FASTA    : Fasta2NonInterLeavedFasta<CR>
menu Phylo.Convert.Wrap\ FASTA      : WrapFasta<CR>

"" DNA
menu Phylo.DNA.RevComp                  : RC<CR>
menu Phylo.DNA.Insert\ Random\ DNA\ Seq : Randseq

"" Run programs
"" need to capture the output from fasttree in a new buffer instead of having
"" to hardcode the output name
menu Phylo.Run.Fasttree\ (DNA)         : ! fasttree -nt % > FaSt.tre<CR><CR>: tabe FaSt.tre<CR>
menu Phylo.Run.Fasttree\ (AA)          : ! fasttree  % > FaSt.tre<CR><CR>: tabe FaSt.tre<CR>
menu Phylo.Run.PAUP                    : ! paup %<CR>
menu Phylo.Run.MrBayes                 : ! mb312 -i %<CR>
"menu Phylo.Phylo.Run\ Phyml :
menu Phylo.Run.MrAIC\ (24\ nt\ models) : ! mraic.pl %<CR><CR>: tabe *.MrAIC.txt<CR>
menu Phylo.Run.MrAIC\ (56\ nt\ models) : ! mraic.pl -modeltest %<CR><CR>: tabe *.MrAIC.txt<CR>
