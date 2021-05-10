let adasd = 0;
let oldlog = console.log;
console.log = function(text) {
  if (adasd >= 2) {
    return;
  }
  oldlog(text);
  ++adasd;
}


/*
 * Copyright 2010, Google Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following disclaimer
 * in the documentation and/or other materials provided with the
 * distribution.
 *     * Neither the name of Google Inc. nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
  
// Marching cubes in Javascript
//
// Yes, this is madness. But this should test those JS engines!
// Does not do simple optimizations like vertex sharing. Nevertheless,
// performance is quite acceptable on Chrome.
//
// Converted from the standard C implementation that's all over the web.

import { vec3 } from "gl-matrix";
import { Simulation } from "./simulation/simulation.js";

function MarchingCubesEffect(resolution) {
  var ext = gl.getExtension("OES_texture_float");
  if (!ext) {
    alert("this machine or browser does not support OES_texture_float");
    return;
  }

  var arrays = tdl.primitives.createCube(1.0)
  var program = tdl.programs.loadProgramFromScriptTags(
      "ray_vs", "ray_fs")
  var textures = [new tdl.textures.ExternalTexture2D()]
  

  var view = new Float32Array(16)

  var model = new tdl.models.Model(program, arrays, textures);

  var eyePosition = new Float32Array([0, 1, 2.1])
  var target = new Float32Array([0, 0, 0])

  // Size of field. 32 is pushing it in Javascript :)
  var size = resolution
  // Deltas
  var delta = 2.0 / size
  var yd = size
  var zd = size * size
  var size3 = size * size * size
  var max_tex_dim = 16384;
  if (size3 > max_tex_dim * 4) {
    alert("Resolution too high! Something's wrong.");
  }

  //var array = cArray(max_tex_dim * 4);

  //var field = array.data
  var field = new Float32Array(max_tex_dim * 4);

  var m4 = tdl.fast.matrix4

  
  var tex = textures[0].texture;
  var tex_level = 0;
  var tex_width = max_tex_dim;
  var tex_height = 1;

  gl.bindTexture(gl.TEXTURE_2D, tex);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

  // Temp buffers used in polygonize.
  var vlist = new Float32Array(12 * 3)
  var nlist = new Float32Array(12 * 3)

  m4.lookAt(view, eyePosition, target, up);

  function lerp(a,b,t) { return a + (b - a) * t; }
  function VIntX(q,pout,nout,offset,isol,x,y,z,valp1,valp2) {
    var mu = (isol - valp1) / (valp2 - valp1);
    pout[offset + 0] = x + mu * delta;
    pout[offset + 1] = y;
    pout[offset + 2] = z;
    nout[offset + 0] = lerp(normal_cache[q],   normal_cache[q+3], mu)
    nout[offset + 1] = lerp(normal_cache[q+1], normal_cache[q+4], mu)
    nout[offset + 2] = lerp(normal_cache[q+2], normal_cache[q+5], mu)
  }
  function VIntY(q,pout,nout,offset,isol,x,y,z,valp1,valp2) {
    var mu = (isol - valp1) / (valp2 - valp1);
    pout[offset + 0] = x;
    pout[offset + 1] = y + mu * delta;
    pout[offset + 2] = z;
    var q2 = q + yd*3
    nout[offset + 0] = lerp(normal_cache[q],   normal_cache[q2], mu)
    nout[offset + 1] = lerp(normal_cache[q+1], normal_cache[q2+1], mu)
    nout[offset + 2] = lerp(normal_cache[q+2], normal_cache[q2+2], mu)
  }
  function VIntZ(q,pout,nout,offset,isol,x,y,z,valp1,valp2) {
    var mu = (isol - valp1) / (valp2 - valp1);
    pout[offset + 0] = x;
    pout[offset + 1] = y;
    pout[offset + 2] = z + mu * delta
    var q2 = q + zd*3
    nout[offset + 0] = lerp(normal_cache[q],   normal_cache[q2], mu)
    nout[offset + 1] = lerp(normal_cache[q+1], normal_cache[q2+1], mu)
    nout[offset + 2] = lerp(normal_cache[q+2], normal_cache[q2+2], mu)
  }

  function compNorm(q) {
    if (normal_cache[q*3] == 0.0) {
      normal_cache[q*3    ] = field[q-1]  - field[q+1]
      normal_cache[q*3 + 1] = field[q-yd] - field[q+yd]
      normal_cache[q*3 + 2] = field[q-zd] - field[q+zd]
    }
  }

  // Returns total number of triangles. Fills triangles.
  // TODO: Optimize to death, add normal calculations so that we can run
  // proper lighting shaders on the results. The grid parameter should be
  // implicit and removed.
  function polygonize(fx, fy, fz, q, isol) {
    var cubeindex = 0;
    var field0 = field[q]
    var field1 = field[q+1]
    var field2 = field[q+yd]
    var field3 = field[q+1+yd]
    var field4 = field[q+zd]
    var field5 = field[q+1+zd]
    var field6 = field[q+yd+zd]
    var field7 = field[q+1+yd+zd]

    if (field0 < isol) cubeindex |= 1;
    if (field1 < isol) cubeindex |= 2;
    if (field2 < isol) cubeindex |= 8;
    if (field3 < isol) cubeindex |= 4;
    if (field4 < isol) cubeindex |= 16;
    if (field5 < isol) cubeindex |= 32;
    if (field6 < isol) cubeindex |= 128;
    if (field7 < isol) cubeindex |= 64;

    // If cube is entirely in/out of the surface - bail, nothing to draw.
    var bits = edgeTable[cubeindex]
    if (bits == 0) return 0;

    var d = delta
    var fx2 = fx + d, fy2 = fy + d, fz2 = fz + d

    // Top of the cube
    if (bits & 1)    {compNorm(q);       compNorm(q+1);       VIntX(q*3,      vlist, nlist, 0, isol, fx,  fy,  fz, field0, field1); }
    if (bits & 2)    {compNorm(q+1);     compNorm(q+1+yd);    VIntY((q+1)*3,  vlist, nlist, 3, isol, fx2, fy,  fz, field1, field3); }
    if (bits & 4)    {compNorm(q+yd);    compNorm(q+1+yd);    VIntX((q+yd)*3, vlist, nlist, 6, isol, fx,  fy2, fz, field2, field3); }
    if (bits & 8)    {compNorm(q);       compNorm(q+yd);      VIntY(q*3,      vlist, nlist, 9, isol, fx,  fy,  fz, field0, field2); }
    // Bottom of the cube
    if (bits & 16)   {compNorm(q+zd);    compNorm(q+1+zd);    VIntX((q+zd)*3,    vlist, nlist, 12, isol, fx,  fy,  fz2, field4, field5); }
    if (bits & 32)   {compNorm(q+1+zd);  compNorm(q+1+yd+zd); VIntY((q+1+zd)*3,  vlist, nlist, 15, isol, fx2, fy,  fz2, field5, field7); }
    if (bits & 64)   {compNorm(q+yd+zd); compNorm(q+1+yd+zd); VIntX((q+yd+zd)*3, vlist, nlist, 18, isol, fx,  fy2, fz2, field6, field7); }
    if (bits & 128)  {compNorm(q+zd);    compNorm(q+yd+zd);   VIntY((q+zd)*3,    vlist, nlist, 21, isol, fx,  fy,  fz2, field4, field6); }
    // Vertical lines of the cube
    if (bits & 256)  {compNorm(q);       compNorm(q+zd);      VIntZ(q*3,        vlist, nlist, 24, isol, fx,  fy,  fz, field0, field4); }
    if (bits & 512)  {compNorm(q+1);     compNorm(q+1+zd);    VIntZ((q+1)*3,    vlist, nlist, 27, isol, fx2, fy,  fz, field1, field5); }
    if (bits & 1024) {compNorm(q+1+yd);  compNorm(q+1+yd+zd); VIntZ((q+1+yd)*3, vlist, nlist, 30, isol, fx2, fy2, fz, field3, field7); }
    if (bits & 2048) {compNorm(q+yd);    compNorm(q+yd+zd);   VIntZ((q+yd)*3,   vlist, nlist, 33, isol, fx,  fy2, fz, field2, field6); }

    cubeindex <<= 4  // Re-purpose cubeindex into an offset into triTable.
    var numtris = 0, i = 0;
    while (triTable[cubeindex + i] != -1) {
      imm.posnormtriv(vlist, nlist,
                      3 * triTable[cubeindex + i + 0],
                      3 * triTable[cubeindex + i + 1],
                      3 * triTable[cubeindex + i + 2])
      i += 3;
      numtris++;
    }
    return numtris;
  }

  function smin(a, b, k) {
    h = Math.min(1, Math.max(0, 0.5+0.5*(b-a)/k));
    return b * (1 - h) + a * h - k*h*(1.0-h);
  }

  function sphereSDF(x, y, z, r) {
    return Math.sqrt(x*x + y*y + z*z) - r;
  }

  // Adds a reciprocal ball (nice and blobby) that, to be fast, fades to zero after
  // a fixed distance, determined by strength and subtract.
  function addBallOld(ballx, bally, ballz, strength, subtract) {
    // Let's solve the equation to find the radius:
    // 1.0 / (0.000001 + radius^2) * strength - subtract = 0
    // strength / (radius^2) = subtract
    // strength = subtract * radius^2
    // radius^2 = strength / subtract
    // radius = sqrt(strength / subtract)
    var radius = size * Math.sqrt(strength / subtract)
    var min_z = Math.floor(ballz * size - radius); if (min_z < 1) {min_z = 1;}
    var max_z = Math.floor(ballz * size + radius); if (max_z > size - 1) max_z = size - 1
    var min_y = Math.floor(bally * size - radius); if (min_y < 1) min_y = 1
    var max_y = Math.floor(bally * size + radius); if (max_y > size - 1) max_y = size - 1
    var min_x = Math.floor(ballx * size - radius); if (min_x < 1) min_x = 1
    var max_x = Math.floor(ballx * size + radius); if (max_x > size - 1) max_x = size - 1
    // Don't polygonize in the outer layer because normals aren't
    // well-defined there.
    for (var z = min_z; z < max_z; z++) {
      var z_offset = size * size * z;
      var fz = z / size - ballz
      var fz2 = fz * fz
      for (var y = min_y; y < max_y; y++) {
        var y_offset = z_offset + size * y;
        var fy = y / size - bally
        var fy2 = fy * fy
        for (var x = min_x; x < max_x; x++) {
          var fx = x / size - ballx
          var val = strength / (0.000001 + fx*fx + fy2 + fz2) - subtract
          if (val > 0.0) field[y_offset + x] += val
        }
      }
    }
  }

  function addBall(x_c, y_c, z_c, radius) {
    for (var z = 0; z < size; ++z) {
      var z_offset = size * size * z;
      var z_w = z / size;

      for (var y = 0; y < size; ++y) {
        var y_offset = z_offset + size * y;
        var y_w = y / size;

        for (var x = 0; x < size; ++x) {
          var x_w = x / size;
          
          field[y_offset + x] = Math.min(field[y_offset + x], sphereSDF(x_w - x_c, y_w - y_c, z_w - z_c, radius), 0.05);
          //if (z == 0)
          //console.log(sphereSDF(x_w - x, y_w - y, z_w - z))
        }
      }
    }
  }

  function addFloor(strength, subtract) {
    var dist = size * Math.sqrt(strength / subtract)
    if (dist > size) dist = size
    for (var y = 0; y < dist; y++) {
      var yy = (y / size) * (y / size)
      var val = strength / (0.0001 + yy) - subtract
      if (val > 0.0) {
        for (var x = 0; x < size; x++)
          for (var z = 0; z < size; z++)
            field[zd * z + y * yd + x] += val
      }
    }
  }

  var firstDraw = true
  var startTime = Date.now() / 1000;
  var lastTime = startTime;

  const gpu = new GPU();
  const sim = new Simulation(gpu, {
    particleDensity: 2000,
    particleBounds: {
      min: vec3.fromValues(0.3, 0.3, 0.3),
      max: vec3.fromValues(0.7, 0.7, 0.7),
    },
    gridBounds: {
      min: vec3.fromValues(0.1, 0.1, 0.1),
      max: vec3.fromValues(0.9, 0.9, 0.9),
    },
  });

  const fillField = gpu.createKernel(function(balls, n, size, radius) {
      let z = Math.floor(this.thread.x / (size * size));
      let y = Math.floor(this.thread.x / size) % size;
      let x = this.thread.x % size;
      var z_w = z / size;
      var y_w = y / size;
      var x_w = x / size;
      let closest = 0;
      let best = 100000;
      // if n too big will need to change loopmaxiterations
      for (let i = 0; i < n; ++i) {
        let cur = (x_w - balls[i][0]) * (x_w - balls[i][0]) + (y_w - balls[i][1]) * (y_w - balls[i][1]) + (z_w - balls[i][2]) * (z_w - balls[i][2]);
        if (cur < best) {
          best = cur;
          closest = i;
        }
      }
      
      //field[y_offset + x] = sphereSDF(x_w - balls[closest][0], y_w - balls[closest][1], z_w - balls[closest][2], radius);
      return Math.sqrt((x_w - balls[closest][0]) * (x_w - balls[closest][0]) + 
             (y_w - balls[closest][1]) * (y_w - balls[closest][1]) + 
             (z_w - balls[closest][2]) * (z_w - balls[closest][2])) - radius;
  }).setLoopMaxIterations(10000).setPipeline(true).setOutput([max_tex_dim * 4]);

  const smooth = gpu.createKernel(function(field, size) {
    let z_c = Math.floor(this.thread.x / (size * size));
    let y_c = Math.floor(this.thread.x / size) % size;
    let x_c = this.thread.x % size;

    let sum = 0;
    let count = 0.001; // weight must be nonzero
    let r = 2;
    for (let z_o = -r; z_o <= r; ++z_o) {
      for (let y_o = -r; y_o <= r; ++y_o) {
        for (let x_o = -r; x_o <= r; ++x_o) {
          let x = x_c + x_o;
          let y = y_c + y_o;
          let z = z_c + z_o;
          
          let w = Math.pow(2.71, -1.5 * Math.sqrt(x_o * x_o + y_o * y_o + z_o * z_o));

          if (x < 0 || x > size-1 || y < 0 || y > size-1 || z < 0 || z > size-1) {
            // cheaper than continue, will try to read invalid data
            w = 0.0;
          }
          
          //let w = x_o == 0 && y_o == 0 && z_o == 0 ? 1.0 : 0.0;
          
          sum += field[z * size * size + y * size + x] * w;
          count += w;
        }
      }
    }

    return sum / count;
  }).setPipeline(true).setOutput([max_tex_dim * 4]);

  this.render = function(framebuffer, time, numblobs) {
    /*m4.perspective(proj, tdl.math.degToRad(60), aspect, 0.1, 500);
    m4.rotationY(world, 0);//time * 0.5)
    m4.translate(world, [0, 0, 0])
    m4.mul(viewproj, view, proj)
    m4.mul(worldview, world, view)
    m4.mul(worldviewproj, world, viewproj)*/

    gl.clearColor(0.0,0.0,0.0,1)
    gl.clearDepth(1.0)
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
    gl.enable(gl.DEPTH_TEST)
    gl.enable(gl.CULL_FACE)

    let curTime = Date.now() / 1000;
    let deltaTime = curTime - lastTime;
    lastTime = curTime;

    let localTime = Date.now() / 1000 - startTime;

    if (deltaTime > 1/20) {
      deltaTime = 1/20;
    }

    //console.log(sim.particles.get(0));
    //console.log(sim.grid.velocityX[0][0][0]);
    sim.step(deltaTime);
    

    var uniformsConst = {
      //u_worldviewproj: worldviewproj,
      //u_worldview: worldview,
      //u_world: world,
      //u_lightDir: [-1.0, 1.0, 1.0],
      //u_lightColor: [116 / 255.0, 204 / 255.0, 244 / 255.0, 1.0],
      //u_ambientUp: [35 / 255.0 * 0.5, 137 / 255.0 * 0.5, 218 / 255.0 * 0.5, 1.0],
      //u_ambientDown: [15 / 255.0 * 0.5, 94 / 255.0 * 0.5, 156 / 255.0 * 0.5, 1.0],
      u_field: textures[0],
      time: localTime
    }

    model.drawPrep(uniformsConst)

    if (firstDraw) {
      // Uncomment to check the speed impact of the field filling.
      //firstDraw = false
      //for (var i = 0; i < size * size * size; i++) {
      // field[i] = 100000.0
      //}
     /*
      let balls = [];
      let n = 30;
      let radius = 0.04;
      
      for (let x = 0; x < n; ++x) {
        for (let z = 0; z < n; ++z) {
          let xp = (x+0.5) / n;
          let zp = (z+0.5) / n;
          let r = Math.sqrt((xp-0.5) * (xp-0.5) + (zp-0.5) * (zp-0.5));
          //let y = 0.1*((Math.sin(40 * r - 1.5*time) + 1) / 2) / Math.abs(10*(Math.max(r, 0.013))) + 0.05;
          let y = 0.3 * Math.pow(Math.cos(10 * r - 1 * time), 2) / Math.max(10*r, 0.5) + 0.05;
          balls.push([xp, y, zp]);
        }
      }*/
      
      let balls = [];
      let radius = 0.05;
      for (let i = 0; i < sim.particles.particleBuffer.length; i += 6) {
        balls.push([sim.particles.particleBuffer[i], sim.particles.particleBuffer[i+1], sim.particles.particleBuffer[i+2]]);
      }
      field = fillField(balls, balls.length, size, radius).toArray()
      //field = smooth(fillField(balls, balls.length, size, radius), size).toArray();
    }

    var isol = 80.0
    imm.begin(gl.TRIANGLES, program);

    /*field[0] = 0.0;
    field[1] = 1.0;
    field[2] = 0.0;
    //field[4] = 0.0;
    field[size + 0] = 0.0;
    field[size + 1] = 0.0;
    field[size + 2] = 1.0;
    field[size * size + 0] = 1.0;
    field[size * size + 1] = 0.0;
    field[size * size + 2] = 0.0;
    field[(size - 1) * size * size + 0] = 0.0;
    field[(size - 1) * size * size + 1] = 1.0;
    field[(size - 1) * size * size + 2] = 0.0;
    
    field[(size - 1) * size * size + (size-1) * size + 0] = 0.0;
    field[(size - 1) * size * size + (size-1) * size + 1] = 0.0;
    field[(size - 1) * size * size + (size-1) * size + 2] = 1.0;

    
    field[(size - 1) * size * size + (size-1) * size + (size-3)] = 1.0;
    field[(size - 1) * size * size + (size-1) * size + (size-2)] = 0.0;
    field[(size - 1) * size * size + (size-1) * size + (size-1)] = 0.0;*/

    /*for (let i = 0; i < size; ++i) {
      for (let j = 0; j < size; ++j) {
        field[i * size + j] = (1 + Math.sin(5*(i/size + j/size))) / 2;
      }
    }*/
    
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texImage2D(gl.TEXTURE_2D, tex_level, gl.RGBA, tex_width, tex_height, 0, gl.RGBA, gl.FLOAT, field);

    gl.activeTexture(gl.TEXTURE0);

    imm.quad2d(-1, -1, 2, 2, 1);

    /*// Triangulate. Yeah, this is slow.
    var size2 = size / 2.0
    for (var z = 2; z < size - 2; z++) {
      var z_offset = size * size * z;
      var fz = (z - size2) / size2 //+ 1
      for (var y = 2; y < size - 2; y++) {
        var y_offset = z_offset + size * y;
        var fy = (y - size2) / size2 //+ 1
        for (var x = 2; x < size - 2; x++) {
          var fx = (x - size2) / size2 //+ 1
          var q = y_offset + x
          polygonize(fx, fy, fz, q, isol)
        }
      }
    }*/
    imm.end();
  }
}

window.MarchingCubesEffect = MarchingCubesEffect;