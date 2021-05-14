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

import { vec3 } from "gl-matrix";
import { Simulation } from "./simulation/simulation.js";

function RayMarchingEffect(resolution, density) {
  var ext = gl.getExtension("OES_texture_float");
  if (!ext) {
    alert("this machine or browser does not support OES_texture_float");
    return;
  }

  var arrays = tdl.primitives.createCube(1.0);
  var program = tdl.programs.loadProgramFromScriptTags("ray_vs", "ray_fs");
  var textures = [new tdl.textures.ExternalTexture2D()];

  var model = new tdl.models.Model(program, arrays, textures);

  var size = resolution;

  var size3 = size * size * size;
  var max_tex_dim = 16384;
  if (size3 > max_tex_dim * 4) {
    alert("Resolution too high! Something's wrong.");
  }

  var field = new Float32Array(max_tex_dim * 4);

  var tex = textures[0].texture;
  var tex_level = 0;
  var tex_width = max_tex_dim;
  var tex_height = 1;

  gl.bindTexture(gl.TEXTURE_2D, tex);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

  var firstDraw = true;
  var startTime = Date.now() / 1000;
  var lastTime = startTime;

  const gpu = new GPU();
  const sim = new Simulation(gpu, {
    particleDensity: density,
    particleBounds: {
      min: vec3.fromValues(0.15, 0.15, 0.15),
      max: vec3.fromValues(0.8, 0.75, 0.35),
    },
    gridBounds: {
      min: vec3.fromValues(0.1, 0.0, 0.1),
      max: vec3.fromValues(0.9, 0.9, 0.9),
    },
  });

  const fillField = gpu
    .createKernel(function (balls, n, size, radius) {
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
        let cur =
          (x_w - balls[i][0]) * (x_w - balls[i][0]) +
          (y_w - balls[i][1]) * (y_w - balls[i][1]) +
          (z_w - balls[i][2]) * (z_w - balls[i][2]);
        if (cur < best) {
          best = cur;
          closest = i;
        }
      }

      return (
        Math.sqrt(
          (x_w - balls[closest][0]) * (x_w - balls[closest][0]) +
            (y_w - balls[closest][1]) * (y_w - balls[closest][1]) +
            (z_w - balls[closest][2]) * (z_w - balls[closest][2])
        ) - radius
      );
    })
    .setLoopMaxIterations(10000)
    .setPipeline(true)
    .setOutput([max_tex_dim * 4]);

  const smooth = gpu
    .createKernel(function (field, size, coefficient) {
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

            // Weighted by e^(-r^2 / c)
            let w = Math.pow(
              2.71,
              (-1 * Math.sqrt(x_o * x_o + y_o * y_o + z_o * z_o)) / coefficient
            );

            if (
              x < 0 ||
              x > size - 1 ||
              y < 0 ||
              y > size - 1 ||
              z < 0 ||
              z > size - 1
            ) {
              // cheaper than continue, will try to read invalid data
              w = 0.0;
            }

            sum += field[z * size * size + y * size + x] * w;
            count += w;
          }
        }
      }

      return sum / count;
    })
    .setPipeline(true)
    .setOutput([max_tex_dim * 4]);

  this.render = function () {
    gl.clearColor(0.0, 0.0, 0.0, 1);
    gl.clearDepth(1.0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.enable(gl.DEPTH_TEST);
    gl.enable(gl.CULL_FACE);

    let curTime = Date.now() / 1000;
    let deltaTime = curTime - lastTime;
    lastTime = curTime;

    let localTime = Date.now() / 1000 - startTime;

    // step the simulation forwards
    deltaTime = Math.min(deltaTime, 1 / 60);
    sim.step(deltaTime / 2);

    var uniformsConst = {
      u_field: textures[0],
      time: localTime,
    };

    model.drawPrep(uniformsConst);

    if (firstDraw) {
      // Set firstDraw = false to only draw 1 frame

      // Sine wave water
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
      }
      */

      let balls = [];
      //let radius = 0.04;
      let radius = window.radiusSlider.value / 100;
      for (let i = 0; i < sim.particles.particleBuffer.length; i += 6) {
        balls.push([
          sim.particles.particleBuffer[i],
          sim.particles.particleBuffer[i + 1],
          sim.particles.particleBuffer[i + 2],
        ]);
      }

      // Swap comment to see with / without smoothing
      //field = fillField(balls, balls.length, size, radius).toArray();
      field = smooth(
        fillField(balls, balls.length, size, radius),
        size,
        Math.max(0.001, window.smoothSlider.value)
      ).toArray();
    }

    // Send the field to GPU, issue draw
    imm.begin(gl.TRIANGLES, program);

    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texImage2D(
      gl.TEXTURE_2D,
      tex_level,
      gl.RGBA,
      tex_width,
      tex_height,
      0,
      gl.RGBA,
      gl.FLOAT,
      field
    );

    gl.activeTexture(gl.TEXTURE0);

    imm.quad2d(-1, -1, 2, 2, 1);
    imm.end();
  };
}

window.RayMarchingEffect = RayMarchingEffect;
