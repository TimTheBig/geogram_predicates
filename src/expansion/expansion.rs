use crate::{Sign, geo_sign};
use core::{cmp::Ordering, fmt};
use smallvec::SmallVec;

#[derive(Clone, Debug)]
pub struct Expansion<const N: usize = 6> {
    /// Starts inline then moves to heap allocation once the inline capacity is exceeded
    data: SmallVec<f64, N>,
}

impl<const N: usize> fmt::Display for Expansion<N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Expansion(len: {}) = ", self.len())?;
        if self.data.spilled() {
            write!(f, "Vec [")?;
        } else {
            write!(f, "[")?;
        }

        if self.len() > 1 {
            for x in self.data.iter().take(self.len() - 1) {
                write!(f, "{x}, ")?;
            }
            write!(f, "{}", self.data.last().expect("len is greater then 1"))?;
        } else if self.len() == 1 {
            write!(f, "{}", self.data.first().expect("check above"))?;
        }
        write!(f, "]")
    }
}

impl<const N: usize> PartialEq for Expansion<N> {
    fn eq(&self, other: &Self) -> bool {
        self.equals(other)
    }
}

impl<const N: usize> PartialOrd for Expansion<N> {
    fn partial_cmp(&self, other: &Expansion<N>) -> Option<Ordering> {
        // Always returns Some(Ordering) because compare() gives a total order.
        let est_self = self.estimate();
        let est_rhs = other.estimate();
        est_self.partial_cmp(&est_rhs)
    }
}

impl<const N: usize> core::ops::Index<usize> for Expansion<N> {
    type Output = f64;
    fn index(&self, idx: usize) -> &f64 {
        &self.data[idx]
    }
}

impl<const N: usize> Default for Expansion<N> {
    fn default() -> Self {
        Self::with_capacity(N)
    }
}

impl<const N: usize> Expansion<N> {
    pub(crate) fn new() -> Self {
        Self::with_capacity(N)
    }

    /// Create a new `Expansion` with the given capacity.
    ///
    /// Internally uses a [`heapless::Vec<f64, N>`](`heapless::Vec`).
    ///
    /// ## Parameters
    /// - `capacity`: maximum number of components, if this is less then the inline capacity it will still be 9.
    ///
    /// ## Examples
    /// ```
    /// # use geogram_predicates::Expansion;
    /// let e: Expansion = Expansion::with_capacity(6);
    /// assert_eq!(e.capacity(), 6); // 6 is the default Expansion capacity
    /// assert_eq!(e.length(), 0);
    /// ```
    /// ```should_panic
    /// # use geogram_predicates::Expansion;
    /// let e = Expansion::<9>::with_capacity(10);
    /// assert_eq!(e.capacity(), 10);
    /// assert_eq!(e.length(), 0);
    /// ```
    pub fn with_capacity(capacity: usize) -> Self {
        debug_assert!(N >= capacity);

        // it would be nice if this used capacity
        Self {
            data: SmallVec::<f64, N>::with_capacity(capacity),
        }
    }

    #[inline]
    pub const fn len(&self) -> usize {
        self.data.len()
    }

    /// Returns `true` if the vector is empty
    #[inline]
    pub const fn is_empty(&self) -> bool {
        self.data.len() == 0
    }

    #[inline]
    pub const fn capacity(&self) -> usize {
        self.data.capacity()
    }

    #[inline]
    pub fn data(&self) -> &[f64] {
        &self.data
    }

    #[inline]
    pub(crate) const fn data_mut(&mut self) -> &mut SmallVec<f64, N> {
        &mut self.data
    }

    pub fn get(&self, i: usize) -> f64 {
        debug_assert!(i < self.data.len());
        self.data[i]
    }

    pub fn get_mut(&mut self, i: usize) -> &mut f64 {
        debug_assert!(i < self.data.len());
        &mut self.data[i]
    }

    /// Assign a single double to this expansion.
    ///
    /// After this call, `self.length() == 1` and `self[0] == a`.
    ///
    /// ## Examples
    /// ```
    /// # use geogram_predicates::Expansion;
    /// let mut e = Expansion::<2>::with_capacity(2);
    /// e.assign(3.14);
    /// assert_eq!(e.length(), 1);
    /// assert_eq!(e[0], 3.14);
    /// ```
    pub fn assign(&mut self, a: f64) -> &mut Self {
        self.data.clear();
        self.data.push(a);

        self
    }

    pub(crate) fn assign_abs(&mut self, rhs: &mut Expansion) -> &mut Self {
        self.data.extend_from_slice(rhs.data_mut());
        for i in 0..rhs.len() {
            self.data[i] = rhs.data[i];
        }

        if self.sign() == Sign::Negative {
            self.negate();
        }
        self
    }

    /// Negate every component of the expansion in place.
    ///
    /// ## Examples
    /// ```
    /// # use geogram_predicates::Expansion;
    /// let mut e = Expansion::from(2.0);
    /// e.negate();
    /// assert_eq!(e[0], -2.0);
    /// ```
    pub fn negate(&mut self) -> &mut Self {
        for v in self.data.iter_mut() {
            *v = -*v;
        }
        self
    }

    pub(crate) fn scale_fast(&mut self, s: f64) -> &mut Self {
        for v in &mut self.data {
            *v *= s;
        }
        self
    }

    /// Estimate the value of this expansion by summing all components.
    ///
    /// This gives a quick—and not fully accurate—“approximate” value.
    ///
    /// ## Examples
    /// ```
    /// # use geogram_predicates::expansion;
    /// let mut e = expansion![1.0, 0.0000001];
    /// assert!(e.estimate() > 1.0);
    /// ```
    pub fn estimate(&self) -> f64 {
        self.data.iter().sum()
    }

    pub fn sign(&self) -> Sign {
        if let Some(data_last) = self.data.last() {
            geo_sign(*data_last)
        } else /* empty */ {
            Sign::Zero
        }
    }

    pub(crate) fn equals(&self, rhs: &Expansion<N>) -> bool {
        self.compare(rhs) == Sign::Zero
    }

    /// Compare two expansions by their estimated values.
    ///
    /// Returns:
    /// - `Ordering::Less` if `self.estimate() < other.estimate()`,
    /// - `Ordering::Equal` if they are (approximately) equal,
    /// - `Ordering::Greater` otherwise.
    ///
    /// ## Examples
    /// ```
    /// # use geogram_predicates::Expansion;
    /// let a = Expansion::from(1.0);
    /// let b = Expansion::from(2.0);
    /// assert!(a < b);
    /// ```
    pub(crate) fn compare(&self, rhs: &Expansion<N>) -> Sign {
        let est_self = self.estimate();
        let est_rhs = rhs.estimate();
        geo_sign(est_self - est_rhs)
    }
}

impl From<f64> for Expansion<1> {
    /// Create a length-1 expansion holding exactly `a`.
    fn from(a: f64) -> Self {
        let mut e = Expansion::with_capacity(1);

        e.data.push(a);
        e
    }
}

impl<const N: usize> From<[f64; N]> for Expansion<N> {
    /// Create an expansion from an array of `N` doubles.
    ///
    /// The resulting expansion has length `N` and
    /// components equal to the array elements in order.
    fn from(arr: [f64; N]) -> Self {
        Expansion {
            // SAFE: arr len is == N
            data: unsafe { SmallVec::from_buf_and_len_unchecked(core::mem::MaybeUninit::new(arr), N) },
        }
    }
}

impl<const N: usize> From<&[f64]> for Expansion<N> {
    /// Create an expansion from a slice of doubles.
    ///
    /// The resulting expansion has length `slice.len()` and
    /// components equal to the slice elements in order.
    fn from(slice: &[f64]) -> Self {
        Expansion {
            data: SmallVec::from_slice(slice),
        }
    }
}

impl<const N: usize> Expansion<N> {
    /// Compute capacity needed to multiply two expansions a and b.
    pub(crate) fn product_capacity<const AN: usize, const BN: usize>(a: &Expansion<AN>, b: &Expansion<BN>) -> usize {
        a.len().saturating_mul(b.len()).saturating_mul(2).max(1)
    }

    /// Assign `self` = a + b (expansion sum).
    /// Naively concatenates the two expansions and then calls `optimize`.
    pub(crate) fn assign_sum<const AN: usize, const BN: usize>(
        &mut self,
        a: &Expansion<AN>,
        b: &Expansion<BN>,
    ) -> &mut Self {
        // todo maybe it's faster if `a` and `b` are optimized first if `self.capacity` is less then `a.len + b.len`
        // if self.capacity() < a.len() + b.len() {
        //     a.optimize();
        //     b.optimize();
        // }

        self.data.extend_from_slice(a.data());
        self.data.extend_from_slice(b.data());

        self.optimize();
        self.data.shrink_to_fit();
        self
    }

    /// Assign `self` = a - b (expansion difference).
    pub(crate) fn assign_diff<const AN: usize, const BN: usize>(
        &mut self,
        a: &Expansion<AN>,
        b: &Expansion<BN>,
    ) -> &mut Self {
        // build negated b in a temp
        let mut nb = Expansion {
            data: b.data.clone(),
        };
        nb.negate();

        self.assign_sum(a, &nb)
    }

    /// Assign `self` = a * b (expansion product).\
    /// Uses Shewchuk's `two_product` on each coefficient pair.
    pub(crate) fn assign_product<const AN: usize, const BN: usize>(
        &mut self,
        a: &Expansion<AN>,
        b: &Expansion<BN>,
    ) -> &mut Self {
        const { assert!(N > 1, "N must be greater then 1") };

        let needed_cap = Self::product_capacity(a, b);

        debug_assert_ne!(self.data.capacity(), 0);
        self.data.resize(needed_cap, 0.0);

        let mut idx = 0;
        // for each pair (i,j), compute two_product and write low,high
        for &ai in &a.data {
            for &bi in &b.data {
                let (low, high) = two_product(ai, bi);
                self.data[idx] = low;
                self.data[idx + 1] = high;
                idx += 2;
            }
        }

        self.optimize();
        self.data.shrink_to_fit();
        self
    }

    /// Remove trailing zero components to maintain a canonical form.
    ///
    /// After this call, `length()` is the smallest index such that
    /// the last component is non-zero, or zero if all components are zero.
    pub(crate) fn optimize(&mut self) -> &mut Self {
        while let Some(&last) = self.data.last() {
            if last == 0.0 {
                self.data.pop();
            } else {
                break;
            }
        }

        self
    }

    /// Compute the capacity needed to form the 3×3 determinant
    /// of the nine expansions a11…a33.
    ///
    /// Mirrors C++:
    /// ```cpp
    /// index_t c11 = det2x2_capacity(a22,a23,a32,a33);
    /// index_t c12 = det2x2_capacity(a21,a23,a31,a33);
    /// index_t c13 = det2x2_capacity(a21,a22,a31,a32);
    /// return 2 * (a11.len()*c11 + a12.len()*c12 + a13.len()*c13);
    /// ```
    pub(crate) fn det3x3_capacity(
        [a11, a12, a13]: [&Expansion<N>; 3],
        [a21, a22, a23]: [&Expansion<N>; 3],
        [a31, a32, a33]: [&Expansion<N>; 3],
    ) -> usize {
        let c11 = Self::det2x2_capacity(a22, a23, a32, a33);
        let c12 = Self::det2x2_capacity(a21, a23, a31, a33);
        let c13 = Self::det2x2_capacity(a21, a22, a31, a32);
        2 * ((a11.len() * c11) + (a12.len() * c12) + (a13.len() * c13))
    }

    /// Assign to `self` the 3×3 determinant of the nine expansions.
    ///
    /// Computes
    /// ```txt
    /// a11·det2x2(a22,a23,a32,a33)
    /// − a12·det2x2(a21,a23,a31,a33)
    /// + a13·det2x2(a21,a22,a31,a32)
    /// ```
    pub(crate) fn assign_det3x3<const IN_N: usize>(
        &mut self,
        [a11, a12, a13]: [&Expansion<IN_N>; 3],
        [a21, a22, a23]: [&Expansion<IN_N>; 3],
        [a31, a32, a33]: [&Expansion<IN_N>; 3],
    ) -> &mut Expansion<N> {
        // 1) build the three 2×2 minors
        let mut m11: Expansion<4> = Expansion::with_capacity(4);
        m11.assign_det2x2(a22, a23, a32, a33);

        let mut m12: Expansion<4> = Expansion::with_capacity(4);
        m12.assign_det2x2(a21, a23, a31, a33);

        let mut m13: Expansion<4> = Expansion::with_capacity(4);
        m13.assign_det2x2(a21, a22, a31, a32);

        // 2) form the three products
        let mut t1: Expansion<4> = Expansion::with_capacity(4);
        t1.assign_product(a11, &m11);

        let mut t2: Expansion<4> = Expansion::with_capacity(4);
        t2.assign_product(a12, &m12);

        let mut t3: Expansion<4> = Expansion::with_capacity(4);
        t3.assign_product(a13, &m13);

        // 3) combine: (t1 - t2) + t3
        let mut tmp = Expansion::<N>::with_capacity(self.capacity());
        tmp.assign_sum(&t1, &t3);
        let res = self.assign_diff(&tmp, &t2);

        res
    }

    /// Compute the capacity needed to form the 2×2 determinant
    /// of the four expansions a11, a12, a21, a22.
    pub(crate) fn det2x2_capacity(
        a11: &Expansion<N>, a12: &Expansion<N>,
        a21: &Expansion<N>, a22: &Expansion<N>,
    ) -> usize {
        Self::product_capacity(a11, a22) + Self::product_capacity(a21, a12)
    }

    /// Assign to `self` the 2×2 determinant of the four expansions:
    ///     a11·a22 − a12·a21
    ///
    /// ## Panics
    /// Panics if `self.capacity()` is less than `det2x2_capacity(...)`.
    pub(crate) fn assign_det2x2<const IN_N: usize>(
        &mut self,
        a11: &Expansion<IN_N>, a12: &Expansion<IN_N>,
        a21: &Expansion<IN_N>, a22: &Expansion<IN_N>,
    ) -> &mut Expansion<N> {
        const { assert!(N > 1, "N must be greater then 1"); }

        // build product a11 * a22
        let mut p1: Expansion<N> = Expansion::new();
        p1.assign_product(a11, a22);

        // build product a12 * a21
        let mut p2: Expansion<N> = Expansion::new();
        p2.assign_product(a12, a21);

        // self = p1 - p2
        self.assign_diff(&p1, &p2)
    }
}

/// Two-product: return (low, high) parts of ai*bi
#[inline]
fn two_product(a: f64, b: f64) -> (f64, f64) {
    let x = a * b;

    (x, two_product_tail(a, b, x))
}

#[inline]
const fn two_product_tail(a: f64, b: f64, x: f64) -> f64 {
    let (ahi, alo) = split(a);
    let (bhi, blo) = split(b);

    let err1 = x - (ahi * bhi);
    let err2 = err1 - (alo * bhi);
    let err3 = err2 - (ahi * blo);

    (alo * blo) - err3
}

const SPLITTER: f64 = 134_217_729f64;

#[inline]
const fn split(a: f64) -> (f64, f64) {
    let c = SPLITTER * a;
    let abig = c - a;
    let ahi = c - abig;
    let alo = a - ahi;

    (ahi, alo)
}

impl<const N: usize> core::ops::Add for Expansion<N> {
    // todo use when generic_const_exprs is stabe
    // type Output = Expansion<{N + N}>;
    type Output = Expansion<N>;

    fn add(self, rhs: Self) -> Self::Output {
        let mut prod: Expansion<N> = Expansion::with_capacity(N + N);
        prod.assign_sum(&self, &rhs);
        prod
    }
}

impl<const N: usize> core::ops::Mul for Expansion<N> {
    // todo use when generic_const_exprs is stabe
    // type Output = Expansion<{N.saturating_mul(b.length()).saturating_mul(2)}>;
    type Output = Expansion<N>;

    fn mul(self, rhs: Expansion<N>) -> Self::Output {
        // todo when generic_const_exprs is stabe Expansion<{SN.max(RN)}>
        let mut prod: Expansion<N> = Expansion::with_capacity(Expansion::<N>::product_capacity(&self, &rhs));

        prod.assign_product(&self, &rhs);
        prod
    }
}
