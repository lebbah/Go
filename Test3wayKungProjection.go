/*
 *
 * Approximate front computation based on Random projection, Sorting and Kung algorithm
 * Date: July 27, 2018
 * Authors:
 *	{Christophe.Cerin,Mustapha.Lebbah,Tarek.Menouer}@lipn.univ-paris13.fr
 *
 */

/*
For controling the font size:

export GNUPLOT_DEFAULT_GDFONT="Aller"

Or add the SetSvg function to go/src/github.com/Arafatk/glot/common.go

func (plot *Plot) SetSvg() error {
	return plot.Cmd(fmt.Sprintf("set terminal svg linewidth 1 font 'Arial,22' size 1280,960 ; set logscale x 2 ; set xtics right logscale;"))
}

ATTENTION : common.go from "github.com/Arafatk/glot" needs also to be
adapted to allow svg outputs and linespoints style

*/

package main

import (
	"container/list"
	//"encoding/json"
	"fmt"
	//"github.com/Arafatk/glot"
	"flag"
	"math"
	"math/rand"
	"reflect"
	"strconv"
	"sync"
	"time"
)

type Value float64
type Matrix [][]Value

/*
 * The distribution to use for calculating the random matrix.
 *  Sparse1 is:
 *    sqrt(3)*{-1 with prob(1/6), 0 with prob(2/3), +1 with prob(1/6)}
 *  Sparse2 is:
 *    {-1 with prob(1/2), +1 with prob(1/2)}
 * Gaussian is:
 *     Gaussian
 */
type RM int

const (
	Sparse1  RM = 0
	Sparse2  RM = 1
	Gaussian RM = 2
)

func GeneRandomMatrix(m1 Matrix, rm RM) (m2 Matrix, ok bool, mmax Value) {
	rows := len(m1)
	if rows < 1 {
		return nil, false, math.MinInt32
	}
	rand.Seed(time.Now().UTC().UnixNano())
	Max := -Value(math.MaxFloat64)
	m2 = make(Matrix, 1)
	m2[0] = make([]Value, rows)
	for k := 0; k < rows; k++ {
		switch rm {
		case Sparse1:
			n := rand.Intn(6)
			if n < 3 {
				m2[0][k] = 1
			} else if n < 4 {
				m2[0][k] = Value(math.Sqrt(3.0))
			} else {
				m2[0][k] = Value(-1 * math.Sqrt(3.0))
			}
		case Sparse2:
			if Value(rand.NormFloat64()) < 0 {
				m2[0][k] = -1
			} else {
				m2[0][k] = 1
			}
		case Gaussian:
			m2[0][k] = Value(rand.NormFloat64() + 3)
		default:
			m2[0][k] = Value(1)
		}
		if m2[0][k] > Value(Max) {
			Max = m2[0][k]
		}
	}
	return m2, true, Max
}

func Multiply(m1, m2 Matrix) (m3 Matrix, ok bool, mmax int) {
	rows, cols, extra := len(m1), len(m2[0]), len(m2)
	if len(m1[0]) != extra {
		return nil, false, -1
	}
	Max := -Value(math.MaxFloat64)
	m3 = make(Matrix, rows)
	for i := 0; i < rows; i++ {
		m3[i] = make([]Value, cols)
		for j := 0; j < cols; j++ {
			for k := 0; k < extra; k++ {
				m3[i][j] += m1[i][k] * m2[k][j]
			}
			if Max < m3[i][j] {
				Max = m3[i][j]
			}
		}
	}
	return m3, true, int(math.Ceil(float64(Max)))
}

func (m Matrix) String() string {
	ib := Value(1)
	rows := len(m)
	cols := len(m[0])
	out := "["
	for r := 0; r < rows; r++ {
		if r > 0 {
			out += ",\n "
		}
		out += "[ "
		for c := 0; c < cols; c++ {
			if c > 0 {
				out += ", "
			}
			if reflect.TypeOf(ib).Kind() == reflect.Float64 {
				out += fmt.Sprintf("%7.3f", m[r][c])
			} else {
				out += fmt.Sprintf("%d", m[r][c])
			}
		}
		out += " ]"
	}
	out += "]"
	return out
}

func InsertionSort(v []int) {
	for j := 1; j < len(v); j++ {
		// Invariant: v[:j] contains the same elements as
		// the original slice v[:j], but in sorted order.
		key := v[j]
		i := j - 1
		for i >= 0 && v[i] > key {
			v[i+1] = v[i]
			i--
		}
		v[i+1] = key
	}
}

// Partition reorders the elements of v so that:
// - all elements in v[:low] are less than p,
// - all elements in v[low:high] are equal to p,
// - all elements in v[high:] are greater than p.
func Partition(v []int, p int) (low, high int) {
	low, high = 0, len(v)
	for mid := 0; mid < high; {
		// Invariant:
		//  - v[:low] < p
		//  - v[low:mid] = p
		//  - v[mid:high] are unknown
		//  - v[high:] > p
		//
		//         < p       = p        unknown       > p
		//     -----------------------------------------------
		// v: |          |          |a            |           |
		//     -----------------------------------------------
		//                ^          ^             ^
		//               low        mid           high
		switch a := v[mid]; {
		case a < p:
			v[mid] = v[low]
			v[low] = a
			low++
			mid++
		case a == p:
			mid++
		default: // a > p
			v[mid] = v[high-1]
			v[high-1] = a
			high--
		}
	}
	return
}

func Pivot(v []int) int {
	n := len(v)
	return Median(v[rand.Intn(n)],
		v[rand.Intn(n)],
		v[rand.Intn(n)])
}

func Median(a, b, c int) int {
	if a < b {
		switch {
		case b < c:
			return b
		case a < c:
			return c
		default:
			return a
		}
	}
	switch {
	case a < c:
		return a
	case b < c:
		return c
	default:
		return b
	}
}

//  Quicksort sorts the elements of v in ascending order.
func Quicksort(v []int) {
	if len(v) < 20 {
		InsertionSort(v)
		return
	}
	p := Pivot(v)
	low, high := Partition(v, p)
	Quicksort(v[:low])
	Quicksort(v[high:])
}

//
//  Helper sorting functions.
//

// is v < w ?
func less(v int, w int) bool {
	if v == w {
		return false
	} else // optimization when reference equal
	{
		return (v < w)
	}
}

// return the index of the median element among a[i], a[j], and a[k]
func median3(a []int, i int, j int, k int) int {
	if less(a[i], a[j]) {
		if less(a[j], a[k]) {
			return j
		} else if less(a[i], a[k]) {
			return k
		} else {
			return i
		}
	} else if less(a[k], a[j]) {
		return j
	} else if less(a[k], a[i]) {
		return k
	} else {
		return i
	}
}

// does v == w ?
func eq(v int, w int) bool {
	if v == w {
		return true
	} // optimization when reference equal
	return false
}

// exchange a[i] and a[j]
func exch(a []int, i int, j int) {
	swap := a[i]
	a[i] = a[j]
	a[j] = swap
}

func MyinsertionSort(a []int, lo int, hi int) {
	i := lo
	for i = lo; i <= hi; i++ {
		j := i
		for j = i; j > lo && less(a[j], a[j-1]); j-- {
			exch(a, j, j-1)
		}
	}
}

//
//
// 3way Quicksort from Robert Sedgewick
// https://www.cs.princeton.edu/~rs/talks/QuicksortIsOptimal.pdf
//

func _3way(a []int, l int, r int) {
	if r <= l {
		return
	}
	v := a[r]
	i := l - 1
	j := r
	p := l - 1
	q := r
	for {
		i++
		for less(a[i], v) {
			i++
		}
		j--
		for less(v, a[j]) {
			if j == l {
				break
			} else {
				j--
			}
		}
		if i >= j {
			break
		}
		exch(a, i, j)
		if eq(a[i], v) {
			p++
			exch(a, p, i)
		}
		if eq(v, a[j]) {
			q--
			exch(a, q, j)
		}
	}
	exch(a, i, r)
	j = i - 1
	i = i + 1
	for k := l; k < p; k, j = k+1, j-1 {
		exch(a, k, j)
	}
	for k := r - 1; k > q; k, i = k-1, i+1 {
		exch(a, k, i)
	}
	_3way(a, l, j)
	_3way(a, i, r)
}

//
// 3way quicksort for arrays of dictionaries/maps
// The dimension for the sort is given as a parameter (dim)
//

// is v < w ?
func lessD(v map[int]int, w map[int]int, dim int) bool {
	if v[dim] == w[dim] {
		return false
	} else // optimization when reference equal
	{
		return (v[dim] < w[dim])
	}
}

// is v > w ?
func gtD(v map[int]int, w map[int]int, dim int) bool {
	if v[dim] == w[dim] {
		return false
	} else // optimization when reference equal
	{
		return (v[dim] > w[dim])
	}
}

// does v == w ?
func eqD(v map[int]int, w map[int]int, dim int) bool {
	if v[dim] == w[dim] {
		return true
	} // optimization when reference equal
	return false
}

// exchange a[i] and a[j]
func exchD(a []map[int]int, i int, j int, dim int) {
	swap := map[int]int{}
	for key, _ := range a[i] {
		swap[key] = a[i][key]
	}
	for key, _ := range a[i] {
		a[i][key] = a[j][key]
	}
	for key, _ := range a[i] {
		a[j][key] = swap[key]
	}
}

func InsertionSortD(v []map[int]int, dim int) {
	for j := 1; j < len(v); j++ {
		// Invariant: v[:j] contains the same elements as
		// the original slice v[:j], but in sorted order.
		key := v[j]
		i := j - 1
		for i >= 0 && gtD(v[i], key, dim) {
			v[i+1] = v[i]
			i--
		}
		v[i+1] = key
	}
}

func _3wayD(a []map[int]int, l int, r int, dim int) {
	//if r <= l {
	//	return
	//}
	if r-l < 4096 {
		InsertionSortD(a[l:r+1], dim)
		return
	}
	v := make(map[int]int)
	for key, value := range a[r] {
		// fmt.Println("Key:", key, "Value:", value)
		v[key] = value
	}
	i := l - 1
	j := r
	p := l - 1
	q := r
	for {
		i++
		for lessD(a[i], v, dim) {
			i++
		}
		j--
		for lessD(v, a[j], dim) {
			if j == l {
				break
			} else {
				j--
			}
		}
		if i >= j {
			break
		}
		exchD(a, i, j, dim)
		if eqD(a[i], v, dim) {
			p++
			exchD(a, p, i, dim)
		}
		if eqD(v, a[j], dim) {
			q--
			exchD(a, q, j, dim)
		}
	}
	exchD(a, i, r, dim)
	j = i - 1
	i = i + 1
	for k := l; k < p; k, j = k+1, j-1 {
		exchD(a, k, j, dim)
	}
	for k := r - 1; k > q; k, i = k-1, i+1 {
		exchD(a, k, i, dim)
	}
	_3wayD(a, l, j, dim)
	_3wayD(a, i, r, dim)
}

func _3wayD_PAR(a []map[int]int, l int, r int, dim int) {
	if r <= l {
		return
	}
	v := make(map[int]int)
	for key, value := range a[r] {
		// fmt.Println("Key:", key, "Value:", value)
		v[key] = value
	}
	i := l - 1
	j := r
	p := l - 1
	q := r
	for {
		i++
		for lessD(a[i], v, dim) {
			i++
		}
		j--
		for lessD(v, a[j], dim) {
			if j == l {
				break
			} else {
				j--
			}
		}
		if i >= j {
			break
		}
		exchD(a, i, j, dim)
		if eqD(a[i], v, dim) {
			p++
			exchD(a, p, i, dim)
		}
		if eqD(v, a[j], dim) {
			q--
			exchD(a, q, j, dim)
		}
	}
	exchD(a, i, r, dim)
	j = i - 1
	i = i + 1
	for k := l; k < p; k, j = k+1, j-1 {
		exchD(a, k, j, dim)
	}
	for k := r - 1; k > q; k, i = k-1, i+1 {
		exchD(a, k, i, dim)
	}

	var wg sync.WaitGroup
	wg.Add(2)
	go func() { _3wayD_PAR(a, l, j, dim); wg.Done() }()
	go func() { _3wayD_PAR(a, i, r, dim); wg.Done() }()
	wg.Wait()
}

func Index(MyTab []map[int]int, nodeCount int, nbDimension int) int {

	var index int = -1
	var stop bool = false

	for i := 0; i < nodeCount; i++ {

		if MyTab[i][nbDimension-1] != -1 && stop == false {
			index = i
			stop = true
		}
	}

	return index
}

func Size(MyTab []map[int]int, nodeCount int, nbDimension int) int {

	var cmpt int = 0

	for i := 0; i < nodeCount; i++ {
		if MyTab[i][nbDimension-1] != -1 {
			cmpt++
		}
	}

	return cmpt
}

func dominate(x map[int]int, y map[int]int, nbDimension int) bool {

	var cmpt_equal int = 0

	for index := 0; index < nbDimension/2; index++ {
		if x[index*2] == y[index*2] {
			cmpt_equal++
		} else {
			if !test_dominate(x[index*2], y[index*2], x[index*2+1]) {
				return false
			}
		}
	}

	//All x and y values are equal => no dominance
	if cmpt_equal == (nbDimension / 2) {
		return false
	}

	return true
}

func test_dominate(x int, y int, min_max int) bool {

	if min_max == 0 {
		return (x < y)
	}
	if min_max == 1 {
		return (x > y)
	}

	fmt.Println("Bad coefficient for min or max objective")
	return false
}

func KUNG(MyTab []map[int]int, nodeCount int, nbDimension int) *list.List {

	/* This function is referenced in: Lixin Ding, Sanyou Zeng and
	Lishan Kang: A fast algorithm on finding the non-dominated set
	in multi-objective optimization.  * Evolutionary Computation,
	2003, Canberra, ACT, Australia, Australia */

	var cmpt, index int = 0, -1
	var stop, stop1 bool = false, false

	var x = make(map[int]int, nbDimension)
	var copie = make(map[int]int, nbDimension)

	MyTabResult := make([]map[int]int, nodeCount)
	CopyMyTab := make([]map[int]int, nodeCount) //Copie of MyTab
	for i := range MyTab {
		MyTabResult[i] = make(map[int]int)
		CopyMyTab[i] = make(map[int]int)
		for j := 0; j < nbDimension; j++ {
			MyTabResult[i][j] = -1
			CopyMyTab[i][j] = MyTab[i][j]
		}
	}

	for {
		index = Index(CopyMyTab, nodeCount, nbDimension) //The first element in CopyMyTab

		if index == -1 {
			break
		}

		for i := 0; i < nbDimension; i++ { //save the first element of CopyMyTab in X
			x[i] = CopyMyTab[index][i]
		}

		if Size(MyTabResult, nodeCount, nbDimension) == 0 { //If MyTabRest is empty

			for i := 0; i < nbDimension; i++ {
				CopyMyTab[index][i] = -1 //Delete the element index (x) form CopyMyTab
			}

			for i := 0; i < nbDimension; i++ {
				//Add the element x to MyTabResult
				MyTabResult[Size(MyTabResult, nodeCount, nbDimension)][i] = x[i]
			}

		} else {
			cmpt = 0
			stop = false

			for {
				if cmpt >= Size(MyTabResult, nodeCount, nbDimension) || stop == true {
					break
				}

				for i := 0; i < nbDimension; i++ {
					copie[i] = MyTabResult[cmpt][i]
				}

				if dominate(copie, x, nbDimension) == true {
					for i := 0; i < nbDimension; i++ {
						CopyMyTab[index][i] = -1 //Delete X from CopyMyTab
					}
					stop = true
				}
				cmpt++
			}

			if stop == false {
				cmpt = 0
				stop1 = false

				for {
					if cmpt >= Size(MyTabResult, nodeCount, nbDimension) || stop1 == true {
						break
					}

					for i := 0; i < nbDimension; i++ {
						copie[i] = MyTabResult[cmpt][i]
					}

					if dominate(x, copie, nbDimension) == true {

						for i := 0; i < nbDimension; i++ {
							MyTabResult[cmpt][i] = -1 //Delete X from CopyMyTab
						}

						stop1 = true
					}
					cmpt++
				}

			}
			if stop == false && stop1 == false {

				for i := 0; i < nbDimension; i++ {
					//Add  X to MyTabResult
					MyTabResult[Size(MyTabResult, nodeCount, nbDimension)][i] = x[i]
				}
				for i := 0; i < nbDimension; i++ {
					CopyMyTab[index][i] = -1 //Delete X from CopyMyTab
				}
			}
		}
	}

	//fmt.Println("Result ")

	l_result := list.New()
	for i := 0; i < nodeCount; i++ {
		if MyTabResult[i][nbDimension-1] != -1 {
			//fmt.Println("Node Id = ", MyTabResult[i][nbDimension-1])
			l_result.PushBack(MyTabResult[i])
		}
	}
	//fmt.Println("\n")

	return l_result

}

func FrontPAR(l *list.List, threshold int) *list.List {

	/* This function is a PARALLEL implementation of Kung’s Algorithm as
	   it is presented in
	   https://engineering.purdue.edu/~sudhoff/ee630/Lecture09.pdf*/

	if l.Len() == 1 {
		return l
	}

	T := list.New()
	B := list.New()

	e := l.Front()
	for i := 0; i < l.Len()/2; i++ {
		var x map[int]int
		x = e.Value.(map[int]int)
		T.PushBack(x)
		e = e.Next()
	}

	for i := l.Len() / 2; i < l.Len(); i++ {
		var x map[int]int
		x = e.Value.(map[int]int)
		B.PushBack(x)
		e = e.Next()
	}

	var wg sync.WaitGroup
	wg.Add(2)

	go func() {
		defer func() { wg.Done() }()
		T = FrontPAR(T, threshold)
	}()
	go func() {
		defer func() { wg.Done() }()
		B = FrontPAR(B, threshold)
	}()
	wg.Wait()

	var size_b int = B.Len()
	var size_t int = T.Len()
	var cmpt int = 0

	eb := B.Front()

	for i := 0; i < size_b; i++ {

		et := T.Front()
		cmpt = 0
		for j := 0; j < size_t; j++ {
			var xt map[int]int
			var xb map[int]int

			xt = et.Value.(map[int]int)
			xb = eb.Value.(map[int]int)

			if dominate(xt, xb, len(xt)) == false {
				cmpt++
			}

			if cmpt == size_t {
				T.PushBack(xb)
			}

			et = et.Next()
		}
		eb = eb.Next()
	}

	return T
}

func Front(l *list.List) *list.List {

	/* This function is a SEQUENTIAL and RECURSIVE implementation of
	   Kung’s Algorithm as it is presented in
	   https://engineering.purdue.edu/~sudhoff/ee630/Lecture09.pdf*/

	if l.Len() == 1 {
		return l
	}

	T := list.New()
	B := list.New()

	e := l.Front()
	for i := 0; i < l.Len()/2; i++ {
		var x map[int]int
		x = e.Value.(map[int]int)
		T.PushBack(x)
		e = e.Next()
	}

	for i := l.Len() / 2; i < l.Len(); i++ {
		var x map[int]int
		x = e.Value.(map[int]int)
		B.PushBack(x)
		e = e.Next()
	}

	T = Front(T)
	B = Front(B)

	var size_b int = B.Len()
	var size_t int = T.Len()
	var cmpt int = 0

	eb := B.Front()

	for i := 0; i < size_b; i++ {

		et := T.Front()
		cmpt = 0
		for j := 0; j < size_t; j++ {
			var xt map[int]int
			var xb map[int]int

			xt = et.Value.(map[int]int)
			xb = eb.Value.(map[int]int)

			if dominate(xt, xb, len(xt)) == false {
				cmpt++
			}

			if cmpt == size_t {
				T.PushBack(xb)
			}

			et = et.Next()
		}
		eb = eb.Next()
	}

	return T
}

func KUNG_recursive(MyTab []map[int]int, nodeCount int, nbDimension int) *list.List {

	l := list.New()

	for i := 0; i < nodeCount; i++ {
		l.PushBack(MyTab[i])
	}

	l = Front(l)

	return l
}

func KUNG_recursive_PAR(MyTab []map[int]int, nodeCount int, nbDimension int) *list.List {

	l := list.New()

	for i := 0; i < nodeCount; i++ {
		l.PushBack(MyTab[i])
	}

	var threshold int
	l = FrontPAR(l, threshold)

	return l
}

func Maximization_objective(MyTab []map[int]int) []map[int]int {

	if MyTab[0][1] == 1 {

		//If function 1 is maximization
		x := make([]int, 5)
		for i := 0; i < len(MyTab)/2; i++ {
			for j := 0; j < 5; j++ {
				x[j] = MyTab[i][j]
			}

			for j := 0; j < 5; j++ {
				MyTab[i][j] = MyTab[len(MyTab)-1-i][j]
			}

			for j := 0; j < 5; j++ {
				MyTab[len(MyTab)-1-i][j] = x[j]
			}
		}
	}
	return MyTab
}

func present(x int, y int, l *list.List) bool {
	e := l.Front()
	for index1 := 0; index1 < l.Len(); index1++ {
		if x == e.Value.(map[int]int)[0] && y == e.Value.(map[int]int)[2] {
			return true
		}
		e = e.Next()
	}
	return false
}

func JaccardIndex(l_exact []*list.List, l_proj []*list.List, taille int, n int, m int) {
	my_mod := int(taille / n)
	fmt.Println("Size of buffers:", taille, "Modulo:", my_mod, "Min:", m)
	result := make([]map[string]float64, my_mod)
	for i := range result {
		result[i] = map[string]float64{"jaccard": 0.0}
	}
	l_index := list.New()
	l_index1 := list.New()
	for index := 0; index < taille; index++ {
		l_index = l_proj[index]
		l_index1 = l_exact[index]
		e := l_index.Front()
		union := l_index.Len()
		intersect := 0
		for index1 := 0; index1 < l_index.Len(); index1++ {
			//fmt.Println(index,index1)
			//for key, value := range e.Value.(map[int]int) {
			//    fmt.Println("Key:", key, "Value:", value)
			//}
			//fmt.Println("Val0: ",e.Value.(map[int]int)[0])
			//fmt.Println("Val1: ",e.Value.(map[int]int)[1])
			//fmt.Println("Val2: ",e.Value.(map[int]int)[2])
			//fmt.Println("Val3: ",e.Value.(map[int]int)[3])
			//fmt.Println("Val4: ",e.Value.(map[int]int)[4])
			if present(e.Value.(map[int]int)[0], e.Value.(map[int]int)[2], l_index1) {
				intersect += 1
			} else {
				union += 1
			}
			e = e.Next()
		}
		result[index%my_mod]["jaccard"] = (result[index%my_mod]["jaccard"] + float64(intersect)/float64(union)) / 2.0
		//fmt.Println(index,"Union",union,"Interset",intersect,"Jaccard Index",result[index%my_mod]["jaccard"])
	}
	// Print results
	for i := range result {
		fmt.Println("Jaccard Index for size", int(float64(m)*math.Pow(float64(2), float64(i))), ":", result[i]["jaccard"])
	}
}

func distance(x int, y int) float64 {
	// Euclidien distance in 2D space
	s := 0.0
	d := float64(x)
	s += d * d
	d = float64(y)
	s += d * d
	return math.Sqrt(s)
}

func JaccardSimilarityCoeff(l_exact []*list.List, l_proj []*list.List, taille int, n int, m int) {
	// We compute the Jaccard Similarity coefficients. We assume that the 2 lists/buffers
	// have equal sizes.
	my_mod := int(taille / n)
	fmt.Println("Size of buffers:", taille, "Modulo:", my_mod, "Min:", m)
	result := make([]map[string]float64, my_mod)
	for i := range result {
		result[i] = map[string]float64{"min": 0.0, "max": 0.0}
	}
	l_index1 := list.New()
	l_index2 := list.New()
	for index := 0; index < taille; index++ {
		l_index1 = l_proj[index]
		l_index2 = l_exact[index]
		e1 := l_index1.Front()
		e2 := l_index2.Front()
		d1 := distance(e1.Value.(map[int]int)[0], e1.Value.(map[int]int)[2])
		d2 := distance(e2.Value.(map[int]int)[0], e2.Value.(map[int]int)[2])
		if d1 > d2 {
			result[index%my_mod]["min"] = result[index%my_mod]["min"] + d2
			result[index%my_mod]["max"] = result[index%my_mod]["max"] + d1
		} else {
			result[index%my_mod]["min"] = result[index%my_mod]["min"] + d1
			result[index%my_mod]["max"] = result[index%my_mod]["max"] + d2
		}
		e1 = e1.Next()
		e2 = e2.Next()
	}
	// Print results
	for i := range result {
		fmt.Println("Jaccard Similarity Coefficient for size", int(float64(m)*math.Pow(float64(2), float64(i))), ":", float64(result[i]["min"]/result[i]["max"]))
	}
}

func main() {
	// Controlling the command line
	//var flagVar1 = flag.Int("projection", 1234, "In order to use the random projection method")
	//var flagVar2 = flag.Int("default", 4321, "In order to use the default method (exact method)")
	flag.Parse()
	if flag.NArg() != 1 {
		fmt.Println("Bag number of arguments\n\tFlag: projection to use the random projection method")
		fmt.Println("\tFlag: default to use the default method")
		fmt.Println("\tFlag: analysis to compute similarity metrics")
		return
	}
	fmt.Println("Arguments: ", flag.Args())
	if flag.Args()[0] == "projection" {
		fmt.Println("Random projection selected")
	} else {
		if flag.Args()[0] == "default" {
			fmt.Println("Default method selected")
		} else {
			if flag.Args()[0] == "analysis" {
				fmt.Println("We compute the similarity metrics of performance")
			} else {
				fmt.Println("Bag argument\n\tFlag: projection to use the random projection method")
				fmt.Println("\tFlag: default to use the default method")
				fmt.Println("\tFlag: analysis to compute similary metrics")
				return
			}
		}
	}
	//if (*flagVar1 == 1) || (*flagVar2) == 2 { return }

	//
	// Parameters of the experiment
	//
	powof2 := int(math.Pow(2.0, 20.0)) // Max size ;
	MinSize := 4096                    // Min size = 4096 by default
	Dim1 := 10000000                   // random values between 0 and Dim1
	Dim2 := 10000000                   // random values between 0 and Dim2
	const N int = 15                   // Number of runs
	const NN int = int(N) + 1          // Number of runs + 1

	fmt.Println("Begin...")
	rand.Seed(time.Now().UTC().UnixNano())
	var result_3way [NN]map[int]float64
	var result_kung [NN]map[int]float64
	var MyTabResult []map[int]int
	// For storing the similarity metrics
	l_exact := make([]*list.List, N*(int(math.Log2(float64(powof2/MinSize)))+1))
	l_proj := make([]*list.List, N*(int(math.Log2(float64(powof2/MinSize)))+1))
	My_ind := 0

	for j := 0; j < int(N)+1; j++ {
		result_3way[j] = make(map[int]float64)
		result_kung[j] = make(map[int]float64)
	}

	for j := 0; j < int(N); j++ {
		for taille := MinSize; taille <= powof2; taille = taille * 2 {
			MyTab := make([]map[int]int, taille)
			l := list.New()
			//l1 := list.New()
			for i := range MyTab {
				MyTab[i] = make(map[int]int)
				MyTab[i][0] = rand.Intn(Dim1) //Value in x for function 1
				MyTab[i][1] = 0               // O if minimization - 1 if maximization
				MyTab[i][2] = rand.Intn(Dim2) //value in y for function 2
				MyTab[i][3] = 0               // O if minimization - 1 if maximization
				MyTab[i][4] = i               //Id value
			}

			fmt.Println("Problem size ", taille)

			dimension := 0
			start := time.Now()
			if flag.Args()[0] == "default" || flag.Args()[0] == "analysis" {
				_3wayD_PAR(MyTab, 0, len(MyTab)-1, dimension)
			}
			if flag.Args()[0] == "projection" || flag.Args()[0] == "analysis" {
				// We transform the input into a matrix
				MatMytab := make(Matrix, 2) // 2 dimensions
				for i := 0; i < 2; i++ {
					MatMytab[i] = make([]Value, taille) // dimension
					for jj := 0; jj < taille; jj++ {
						MatMytab[i][jj] = Value(MyTab[jj][2*i])
						//MatMytab[i][jj] = Value(MyTab[jj][2])
					}
				}
				//fmt.Printf("Initial input (under matrix form):\n%s\n\n", MatMytab)

				// We generare the Random Matrix for the projection
				DD, ok, Max := GeneRandomMatrix(MatMytab, Gaussian)
				//fmt.Printf("max = %7.3f\n", Max)
				if !ok && Max < 1 {
					panic("Invalid dimensions")
				}

				// We compute the product i.e. the projection
				PPPP, ok, mm := Multiply(DD, MatMytab)
				if !ok {
					panic("Invalid dimensions DD x MatMytab")
				}
				//fmt.Printf("Projection -- Product of DD and MatMytab (%d):\n%s\n\n", mm,PPPP)

				// We reorganize the elements based on information coming from PPPP
				MyTabInter := make([]int, mm+1, mm+1)
				for i := range MyTabInter {
					MyTabInter[i] = -1
				}
				NB := 0
				jj := 0
				for i := range PPPP[0] {
					jj = int(math.Ceil(float64(PPPP[0][i])))
					if jj > mm || jj < 0 {
						fmt.Println(len(PPPP[0]), " ", i, " ", jj, " ", mm)
						fmt.Println("A coefficient of the Random projection matrix is negative => change your Gaussian function")
						if jj < 0 {
							jj = -jj
						}
						if jj > mm {
							jj = mm
						}
					}
					if MyTabInter[jj] == -1 {
						NB++
					}
					MyTabInter[jj] = i
				}
				//fmt.Printf("NB elements %d over %d\n",NB,mm+1)

				// We build the final result
				MyTabResult = make([]map[int]int, NB)
				jj = 0
				for i := range MyTabInter {
					if MyTabInter[i] != -1 {
						MyTabResult[jj] = make(map[int]int)
						//fmt.Printf("%d %d\n",MyTabInter[i], mm)
						MyTabResult[jj][0] = int(MatMytab[0][MyTabInter[i]]) //Value in x for function 1
						MyTabResult[jj][1] = 0                               // O if minimization - 1 if maximization
						MyTabResult[jj][2] = int(MatMytab[1][MyTabInter[i]]) //value in y for function 2
						MyTabResult[jj][3] = 0                               // O if minimization - 1 if maximization
						MyTabResult[jj][4] = jj                              //Id value
						jj++
					}
				}
				// We sort the projection
				_3wayD_PAR(MyTabResult, 0, len(MyTabResult)-1, dimension)
			}
			elapsed := time.Since(start)
			/*
				// Check if the first component is sorted
				ss := len(MyTabResult) - 1
				for i := range MyTabResult[0:ss]{
				    	if  MyTabResult[i][0] > MyTabResult[i+1][0] {
					   fmt.Println("Not sorted at position: ",i," over ",ss+1," elements"," (",MyTabResult[i][0],",",MyTabResult[i+1][0],")")
					   break
					} else {
					  fmt.Print(MyTabResult[i][0]," ")
					}
				}
				fmt.Println("")
			*/
			if flag.Args()[0] == "default" {
				fmt.Printf("3way Quicksort parallel takes %s\n", elapsed)
			} else {
				fmt.Printf("Random projection takes %s (Number of elements: %d)\n", elapsed, len(MyTabResult))
			}
			result_3way[j][taille] = float64(elapsed / time.Millisecond)
			result_3way[N][taille] = 0.0
			if flag.Args()[0] == "default" {
				MyTab = Maximization_objective(MyTab)
			} else {
				MyTabResult = Maximization_objective(MyTabResult)
			}

			//l = KUNG(MyTab, len(MyTab), 5)
			start = time.Now()
			if flag.Args()[0] == "analysis" {
				l_exact[My_ind] = KUNG_recursive_PAR(MyTab, len(MyTab), 5)
				l_proj[My_ind] = KUNG_recursive_PAR(MyTabResult, len(MyTabResult), 5)
				My_ind++
			} else {
				if flag.Args()[0] == "default" {
					l = KUNG_recursive_PAR(MyTab, len(MyTab), 5)
				} else {
					l = KUNG_recursive_PAR(MyTabResult, len(MyTabResult), 5)
				}
			}
			elapsed = time.Since(start)
			if flag.Args()[0] != "analysis" {
				fmt.Printf("Kung parallel takes %s (Number of elements in the front: %d)\n", elapsed, l.Len())
			}
			result_kung[j][taille] = float64(elapsed / time.Millisecond)
			result_kung[N][taille] = 0.0
		}
	}

	// We print a summary of results
	fmt.Println("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-")
	if flag.Args()[0] == "projection" {
		fmt.Println("Mean execution time (N experiments) for Random projection")
	} else {
		if flag.Args()[0] == "default" {
			fmt.Println("Mean execution time (N experiments) for 3way")
		}
	}
	if flag.Args()[0] != "analysis" {
		for j := 0; j < int(N); j++ {
			for key, value := range result_3way[j] {
				//fmt.Println("Key:", key, "Value:", value)
				result_3way[N][int(key)] += float64(value)
			}
		}
		for key, value := range result_3way[N] {
			fmt.Println("Key:", key, "Value:", value/float64(N))
		}

		fmt.Println("Mean execution time (N experiments) for Kung")
		for j := 0; j < int(N); j++ {
			//s := 0.0
			for key, value := range result_kung[j] {
				//fmt.Println("Key:", key, "Value:", value)
				result_kung[N][int(key)] += float64(value)
			}
		}
		for key, value := range result_kung[N] {
			fmt.Println("Key:", key, "Value:", value/float64(N))
		}
	}
	/*
		dimensions := 2
		// The dimensions supported by the plot
		persist := false
		debug := false
		plot, _ := glot.NewPlot(dimensions, persist, debug)
		pointGroupName := "Kung"
	*/
	points := [][]float64{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}
	points_orig := [][]float64{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}

	/*
		points := [][]float64{{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}}
		points_orig := [][]float64{{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}}
	*/

	i := 0
	max1 := 0.0
	for key := MinSize; key <= powof2; key *= 2 {
		points[0][i] = float64(key)
		points[1][i] = result_3way[N][key] / float64(N)
		i++
		if result_3way[N][key]/float64(N) > max1 {
			max1 = result_3way[N][key] / float64(N)
		}
	}

	i = 0
	max2 := 0.0
	for key := MinSize; key <= powof2; key *= 2 {
		points_orig[0][i] = float64(key)
		points_orig[1][i] = result_kung[N][key] / float64(N)
		i++
		if result_kung[N][key]/float64(N) > max2 {
			max2 = result_kung[N][key] / float64(N)
		}
	}

	var maxmax = 0.0
	if max2 > max1 {
		maxmax = max2
	} else {
		maxmax = max1
	}

	/*
		//fmt.Println(maxmax)
		//fmt.Println(points_orig)
		//plot.SetLogscale("x", 2)
		// Adding a point group
		plot.AddPointGroup(pointGroupName, "linespoints", points_orig)
		plot.AddPointGroup("3way", "linespoints", points)
		// A plot type used to make points/ curves and customize and save them as an image.
		plot.SetTitle("Performance of 3way and Kung")
		// Optional: Setting the title of the plot
		plot.SetXLabel("Input size (tuples)")
		plot.SetYLabel("Time (ms)")
		// Optional: Setting label for X and Y axis
		//plot.SetXrange(MinSize, powof2)
		plot.SetYrange(0, int(maxmax*1.115))
		// Optional: Setting axis ranges
		// Bug on MAcOs with SetFormat() => doesn't work because pdf is obsolet
		// Gnuplot uses pdfcairo; and because svg is not recognized with gplot
		// We need to customize common.go of gplot distribution
		plot.SetFormat("svg")
		//plot.SetSvg()
		var str string
		str = "./front" + strconv.FormatFloat(float64(time.Now().UTC().UnixNano()), 'f', -1, 64) + ".svg"
		fmt.Println("Picture name: ", str)
		plot.SavePlot(str)
	*/
	if flag.Args()[0] != "analysis" {
		fmt.Println("Print on the screen a Gnuplot program to execute on the command line")
		str := "\"front" + strconv.FormatFloat(float64(time.Now().UTC().UnixNano()), 'f', -1, 64) + ".svg\""
		fmt.Println("Picture name: ", str)
		fmt.Println("set terminal svg linewidth 1 font 'Arial,30' size 1280,960 ;")
		if flag.Args()[0] == "projection" {
			fmt.Println("set title \"Kung method = 3way + Random projection\"")
		} else {
			fmt.Println("set title \"Kung method = 3way + front\"")
		}
		fmt.Println("set logscale x 2 ;set format x \"2^{%L}\" ; set xtics center logscale;")
		fmt.Println("set output ", str)
		fmt.Printf("set yrange [0:%d]\n", int(maxmax*1.115))
		fmt.Println("set xlabel \"Input size (tuples)\"")
		fmt.Println("set ylabel \"Time (ms)\"")
		if flag.Args()[0] == "projection" {
			fmt.Println("plot '-' using 1:2 with linespoints title \"Rand. proj\", '-' using 1:2 with linespoints title \"Front\" ")
		} else {
			fmt.Println("plot '-' using 1:2 with linespoints title \"3way\", '-' using 1:2 with linespoints title \"Front\" ")
		}
		for key := MinSize; key <= powof2; key *= 2 {
			fmt.Println(key, " ", result_3way[N][key]/float64(N))
		}
		fmt.Println("e")
		for key := MinSize; key <= powof2; key *= 2 {
			fmt.Println(key, " ", result_kung[N][key]/float64(N))
		}
		fmt.Println("e")
	} else {
		//fmt.Println("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-")
		fmt.Println("We compute now the Jaccard indexes")
		JaccardIndex(l_exact, l_proj, N*(int(math.Log2(float64(powof2)/float64(MinSize)))+1), N, MinSize)
		fmt.Println("We compute now the Jaccard Similarity Coefficients")
		JaccardSimilarityCoeff(l_exact, l_proj, N*(int(math.Log2(float64(powof2)/float64(MinSize)))+1), N, MinSize)
	}

}
