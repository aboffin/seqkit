// Copyright © 2016 Wei Shen <shenwei356@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package cmd

import (
	"fmt"
	"io"
	"runtime"
	"sort"
	"github.com/dustin/go-humanize"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/util/math"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
	"github.com/tatsushid/go-prettytable"
)

// statCmd represents the stat command
var statCmd = &cobra.Command{
	Use:   "stat",
	Short: "simple statistics of FASTA files",
	Long: `simple statistics of FASTA files

`,
	Run: func(cmd *cobra.Command, args []string) {
		config := getConfigs(cmd)
		alphabet := config.Alphabet
		idRegexp := config.IDRegexp
		outFile := config.OutFile
		seq.AlphabetGuessSeqLenghtThreshold = config.AlphabetGuessSeqLength
		seq.ValidateSeq = false
		runtime.GOMAXPROCS(config.Threads)

		files := getFileList(args)

		outfh, err := xopen.Wopen(outFile)
		checkError(err)
		defer outfh.Close()

		var fastxReader *fastx.Reader
		var num, l, lenMin, lenMax, lenSum, lenN50 uint64
		var seqlens []uint64
		var a ctglen
		var seqFormat, t string
		statInfos := []statInfo{}
		for _, file := range files {
			fastxReader, err = fastx.NewReader(alphabet, file, idRegexp)
			checkError(err)

			seqFormat = ""
			num, lenMin, lenMax, lenSum, lenN50 = 0, ^uint64(0), 0, 0, 0
			for {
				record, err := fastxReader.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(err)
					break
				}
				num++
				if seqFormat == "" {
					if len(record.Seq.Qual) > 0 {
						seqFormat = "FASTQ"
					} else {
						seqFormat = "FASTA"
					}
				}
				l = uint64(len(record.Seq.Seq))
				lenSum += l
				if l < lenMin {
					lenMin = l
				}
				if l > lenMax {
					lenMax = l
				}
				seqlens=append(seqlens, l)
				a = seqlens
				sort.Sort(sort.Reverse(a))
			}
			csum := make(ctglen, len(a))
			csum[0] = a[0]	
			i := 1
			for i < len(a) {
				csum[i] = a[i] + csum[i-1]
				i++
			}
			for i, clen := range(csum) {
				fmt.Println("len: ", a[i], ", clen: ", clen, ", lenSum/2: ", lenSum/2)
				if clen>= (lenSum/2) {
					lenN50 = a[i]
					break
				}
			}
			if fastxReader.Alphabet() == seq.DNAredundant {
				t = "DNA"
			} else if fastxReader.Alphabet() == seq.RNAredundant {
				t = "RNA"
			} else {
				t = fmt.Sprintf("%s", fastxReader.Alphabet())
			}

			if num == 0 {
				statInfos = append(statInfos, statInfo{file, seqFormat, t,
					int64(num), int64(lenSum), int64(lenMin),
					0, int64(lenMax), int64(lenN50)})

			} else {
				statInfos = append(statInfos, statInfo{file, seqFormat, t,
					int64(num), int64(lenSum), int64(lenMin),
					math.Round(float64(lenSum)/float64(num), 1), 
					int64(lenMax), int64(lenN50)})
			}
		}

		// format output
		tbl, err := prettytable.NewTable([]prettytable.Column{
			{Header: "file"},
			{Header: "format"},
			{Header: "type"},
			{Header: "num_seqs", AlignRight: true},
			{Header: "sum_len", AlignRight: true},
			{Header: "min_len", AlignRight: true},
			{Header: "avg_len", AlignRight: true},
			{Header: "max_len", AlignRight: true},
			{Header: "n50_len", AlignRight: true}}...)
		checkError(err)
		tbl.Separator = "  "

		for _, info := range statInfos {
			tbl.AddRow(
				info.file,
				info.format,
				info.t,
				humanize.Comma(info.num),
				humanize.Comma(info.lenSum),
				humanize.Comma(info.lenMin),
				humanize.Commaf(info.lenAvg),
				humanize.Comma(info.lenMax),
				humanize.Comma(info.lenN50))
		}
		outfh.Write(tbl.Bytes())
	},
}

type statInfo struct {
	file   string
	format string
	t      string
	num    int64
	lenSum int64
	lenMin int64
	lenAvg float64
	lenMax int64
	lenN50 int64
}
// define sort.interface
type ctglen []uint64
func (a ctglen) Len() int           { return len(a) }
func (a ctglen) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a ctglen) Less(i, j int) bool { return a[i] < a[j] }

func init() {
	RootCmd.AddCommand(statCmd)
}
