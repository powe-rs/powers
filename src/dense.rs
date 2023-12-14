use derive_builder::Builder;
use num_traits::{One, Zero};
use std::ops::{Add, AddAssign, Mul};

#[derive(Default, Builder)]
#[builder(default, build_fn(name = "pre_build", validate = "Self::validate"))]
pub struct Mat<T: Default> {
    #[builder(setter(custom))]
    rows: usize,
    #[builder(setter(custom))]
    cols: usize,

    values: Vec<T>,

    // Matrix element values are stored in row-major order (C-style)
    // or column-major order (Fortran-style).
    #[builder(setter(custom))]
    col_major: bool,
}

impl<T> MatBuilder<T>
where
    T: Default,
{
    pub fn col_major(&mut self) -> &mut Self {
        self.col_major = Some(true);
        self
    }

    pub fn row_major(&mut self) -> &mut Self {
        self.col_major = Some(false);
        self
    }

    pub fn ones(&mut self) -> &mut Self
    where
        T: Clone + One,
    {
        self.values = Some(vec![T::one(); self.rows.unwrap() * self.cols.unwrap()]);
        self
    }

    pub fn build(&self) -> Result<Mat<T>, MatBuilderError>
    where
        T: Clone + Zero,
    {
        let mut mat = self.pre_build()?;
        if self.values.is_none() {
            mat.values = vec![T::zero(); mat.rows * mat.cols];
        }
        Ok(mat)
    }

    fn validate(&self) -> Result<(), String> {
        if let Some(values) = &self.values {
            let expect = self.rows.unwrap() * self.cols.unwrap();
            if values.len() != expect {
                return Err(format!(
                    "values length ({}) must be rows * cols ({} * {} = {})",
                    values.len(),
                    self.rows.unwrap(),
                    self.cols.unwrap(),
                    expect
                ));
            }
        }
        Ok(())
    }
}

impl<T> Mat<T>
where
    T: Default + Copy,
{
    pub fn new(rows: usize, cols: usize) -> MatBuilder<T>
    where
        T: Clone + Zero,
    {
        MatBuilder {
            rows: Some(rows),
            cols: Some(cols),
            // values: Some(vec![T::zero(); rows * cols]),
            ..Default::default()
        }
    }

    pub fn identity(n: usize) -> MatBuilder<T>
    where
        T: Clone + Zero + One,
    {
        let mut values = vec![T::zero(); n * n];
        for i in 0..n {
            values[i * n + i] = T::one();
        }
        MatBuilder {
            rows: Some(n),
            cols: Some(n),
            values: Some(values),
            ..Default::default()
        }
    }

    pub fn with_diagonal(diag: Vec<T>) -> MatBuilder<T>
    where
        T: Clone + Zero + One,
    {
        let n = diag.len();
        let mut values = vec![T::zero(); n * n];
        for i in 0..n {
            values[i * n + i] = diag[i];
        }
        MatBuilder {
            rows: Some(n),
            cols: Some(n),
            values: Some(values),
            ..Default::default()
        }
    }

    pub fn rows(&self) -> usize {
        self.rows
    }
    pub fn cols(&self) -> usize {
        self.cols
    }
    pub fn shape(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }

    pub fn values(&self) -> &Vec<T> {
        &self.values
    }

    pub fn values_mut(&mut self) -> &mut Vec<T> {
        &mut self.values
    }

    #[inline]
    fn ix(&self, row: usize, col: usize) -> usize {
        if !self.col_major {
            row * self.cols + col
        } else {
            col * self.rows + row
        }
    }

    pub fn get_ref(&self, row: usize, col: usize) -> &T {
        assert!(row < self.rows);
        assert!(col < self.cols);
        &self.values[self.ix(row, col)]
    }

    pub fn get_ref_mut(&mut self, row: usize, col: usize) -> &mut T {
        assert!(row < self.rows);
        assert!(col < self.cols);
        let i = self.ix(row, col);
        &mut self.values[i]
    }

    #[inline]
    pub fn get(&self, row: usize, col: usize) -> T
    where
        T: Copy,
    {
        self.values[self.ix(row, col)]
    }

    #[inline]
    pub fn set(&mut self, row: usize, col: usize, v: T) {
        let i = self.ix(row, col);
        self.values[i] = v
    }

    pub fn row(&self, row: usize) -> impl Iterator<Item = &T> {
        assert!(row < self.rows);
        (0..self.cols).map(move |col| self.get_ref(row, col))
    }

    pub fn col(&self, col: usize) -> impl Iterator<Item = &T> {
        assert!(col < self.cols);
        (0..self.rows).map(move |row| self.get_ref(row, col))
    }

    // pub fn row_mut(&mut self, row: usize) -> impl Iterator<Item = &mut T> {
    //     assert!(row < self.rows);
    //     (0..self.cols).map(move |col| self.get_ref_mut(row, col))
    //     // (0..self.cols).map(move |col| {
    //     //     let i = self.ix(row, col);
    //     //     &mut self.values[i]
    //     // })
    // }

    // pub fn col_mut(&mut self, col: usize) -> impl Iterator<Item = &mut T> {
    //     assert!(col < self.cols);
    //     (0..self.rows).map(move |row| self.get_ref_mut(row, col))
    // }

    pub fn select_rows(&self, rows: &[usize]) -> Self
    where
        T: Clone,
    {
        let mut values = Vec::with_capacity(rows.len() * self.cols);
        if self.col_major {
            for &r in rows {
                for c in 0..self.cols {
                    values.push(self.get(r, c));
                }
            }
        } else {
            for &r in rows {
                let i = self.ix(r, 0);
                let j = self.ix(r, self.cols - 1);
                values.extend(self.values[i..=j].iter())
            }
        }
        Self {
            rows: rows.len(),
            cols: self.cols,
            values,
            col_major: self.col_major,
        }
    }

    pub fn select_cols(&self, cols: &[usize]) -> Self
    where
        T: Clone,
    {
        let mut values = Vec::with_capacity(self.rows * cols.len());
        if self.col_major {
            for &c in cols {
                for r in 0..self.rows {
                    values.push(self.get(r, c));
                }
            }
        } else {
            for r in 0..self.rows {
                for &c in cols {
                    values.push(self.get(r, c));
                }
            }
        }
        Self {
            rows: self.rows,
            cols: cols.len(),
            values,
            col_major: self.col_major,
        }
    }

    pub fn diagonal(&self) -> impl Iterator<Item = &T> {
        assert_eq!(self.rows, self.cols);
        (0..self.rows).map(move |i| self.get_ref(i, i))
    }

    pub fn mat_vec(&self, b: &[T]) -> Vec<T>
    where
        T: Mul<Output = T> + Add<Output = T> + Zero + Copy,
    {
        // mat_vec(self.rows, self.cols, &self.data, b, true)

        assert_eq!(b.len(), self.cols);

        let mut y = Vec::with_capacity(b.len());
        for i in 0..self.rows {
            y.push(dot(&self.values[self.ix(i, 0)..], b))
        }
        y
    }

    pub fn mat_mat(&self, b: &Self) -> Self
    where
        T: Mul<Output = T> + AddAssign + Zero + Copy,
    {
        assert_eq!(
            self.cols, b.rows,
            "rows of b {} must equal columns of a {}",
            b.rows, self.cols
        );
        // assert_eq!(self.col_major, false);
        // assert_eq!(b.col_major, false);

        let mut c = vec![T::zero(); self.rows * b.cols];

        for i in 0..self.rows {
            for j in 0..b.cols {
                let c_ij = self.ix(i, j);
                for k in 0..self.cols {
                    let ax_ik = self.get(i, k);
                    let bx_kj = b.get(k, j);
                    c[c_ij] += ax_ik * bx_kj;
                }
            }
        }
        Self {
            rows: self.rows,
            cols: b.cols,
            values: c,
            col_major: self.col_major, // TODO
        }
    }
}

/// Returns a vector with the indexes of the nonzero elements of `a`.
pub fn find(a: &[f64]) -> Vec<usize> {
    a.iter()
        .enumerate()
        .filter_map(|(i, v)| if *v != 0.0 { Some(i) } else { None })
        .collect()
}

/// Computes the dot-product of `a` and `b`.
pub fn dot<T>(a: &[T], b: &[T]) -> T
where
    T: Mul<Output = T> + Add<Output = T> + Zero + Copy,
{
    return a
        .iter()
        .zip(b)
        .map(|(&ai, &bi)| ai * bi)
        .reduce(|x, y| x + y)
        .unwrap_or(T::zero());
}
