/**
 * Create a new 3D array of zeros with the given dimensions.
 *
 * @param {number} x x dimension length
 * @param {number} y y dimension length
 * @param {number} z z dimension length
 * @returns The new array, all values set to zero.
 */
const initialize3DArray = (x, y, z) => {
  let a = [];
  for (let i = 0; i < x; i++) {
    a.push([]);
    for (let j = 0; j < y; j++) {
      a[i].push([]);
      for (let k = 0; k < z; k++) {
        a[i][j].push(0.0);
      }
    }
  }
  return a;
};

export { initialize3DArray };
