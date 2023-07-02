class Node {
    constructor(x, y, neighbors) {
        this.x = x;
        this.y = y;
        this.neighbors = neighbors;
        this.pairs = {};
        this.points = [];
    }
    
    expandPoints(nodes, thickness, pointsSet) {
        const angles = this.getAngles(nodes);
        if (angles.length > 1) {
            angles.push([angles[0][0] + 2 * Math.PI, angles[0][1]]);
            this.pairs[angles[0][1]] = [];
            const temp = pointsSet.length;
            for (let i = 0; i < angles.length - 1; i++) {
                this.pairs[angles[i][1]].push(pointsSet.length);
                this.pairs[angles[i + 1][1]] = [pointsSet.length];
                this.points.push(pointsSet.length);
                pointsSet.push(this.getPoint(angles[i][0], angles[i + 1][0], thickness));
            }
            this.pairs[angles[0][1]].push(temp);
        } else if (angles.length === 1) {
            this.pairs[angles[0][1]] = [pointsSet.length, pointsSet.length + 1];
            this.points.push(pointsSet.length);
            this.points.push(pointsSet.length + 1);
            pointsSet.push(this.getPoint(angles[0][0] - Math.PI, angles[0][0], thickness));
            pointsSet.push(this.getPoint(angles[0][0], angles[0][0] + Math.PI, thickness));
        }
    }
    
    getAngles(nodes) {
    const angles = [];
    for (const n of this.neighbors) {
        let angle = 0.0;
        const dx = nodes[n].x - this.x;
        const dy = nodes[n].y - this.y;
        if (dx === 0) {
            if (dy > 0) {
                angle = Math.PI / 2;
            } else if (dy < 0) {
                angle = 3 * Math.PI / 2;
            } else {
                angle = 0;
            }
        } else if (dx > 0 && dy >= 0) {
            angle = Math.atan(dy / dx);
        } else if (dx < 0 && dy > 0) {
            angle = Math.PI / 2 + Math.atan(-dx / dy);
        } else if (dx < 0 && dy <= 0) {
            angle = Math.PI + Math.atan(dy / dx);
        } else if (dx > 0 && dy < 0) {
            angle = 3 * Math.PI / 2 + Math.atan(-dx / dy);
        }
        angles.push([angle, n]);
    }
    angles.sort((a, b) => a[0] - b[0]);
    return angles;
    }
    
    getPoint(angle1, angle2, thickness) {
        const angle = (angle1 + angle2) / 2;
        const length = thickness / Math.sin((angle2 - angle1) / 2);
        const x = this.x + length * Math.cos(angle);
        const y = this.y + length * Math.sin(angle);
        return [x, y];
    }
}

class Wall {
    constructor(data, base_index, offset, normal_offset) {
    this.thickness = data['thickness'];
    this.height = data['height'];
    this.nodes = {};
    this.points = [];
    this.normals = [];
    this.normal_map = {};
    this.faces = [];
    this.base_index = base_index;
    this.offset = offset;
    this.normal_offset = normal_offset;
    this.loadNodes(data['vertex']);
    }
    loadNodes(vertexs) {
        for (let i = 0; i < vertexs.length; i++) {
            this.nodes[i + this.base_index] = new Node(vertexs[i]['x'], vertexs[i]['y'], vertexs[i]['neighbor']);
        }
    }
    
    reconstruct() {
        this.expandNodes();
        this.rebuild();
        this.NodeProcess();
        return [Object.keys(this.nodes).length, 2 * this.points.length, this.normals.length];
    }
    
    expandNodes() {
        for (const n of Object.values(this.nodes)) {
            n.expandPoints(this.nodes, this.thickness / 2, this.points);
        }
    }
    
    computeNormal(points) {
        const x1 = points[0][0];
        const y1 = points[0][1];
        const x2 = points[1][0];
        const y2 = points[1][1];
        const vector = [x2 - x1, y2 - y1];
        let normal = [-vector[1], vector[0]];
        const magnitude = Math.sqrt(normal[0] ** 2 + normal[1] ** 2);
        normal = [normal[0] / magnitude, normal[1] / magnitude];
        const key = [normal[0].toFixed(4), normal[1].toFixed(4), 0.0];
        if (!this.normal_map[key]) {
            this.normal_map[key] = this.normals.length;
            this.normals.push(key);
        }
        return this.normal_map[key];
    }
    
    rebuild() {
        for (let node in this.nodes) {
          for (let neighbor of this.nodes[node].neighbors) {
            let firstTop = this.nodes[node].pairs[neighbor][0];
            let secondTop = this.nodes[node].pairs[neighbor][1];
            let neighborFirstTop = this.nodes[neighbor].pairs[node][0];
            let neighborSecondTop = this.nodes[neighbor].pairs[node][1];
            if (node < neighbor) {
              if (!([0.0, 0.0, 1.0] in this.normal_map)) {
                this.normal_map[[0.0, 0.0, 1.0]] = this.normals.length;
                this.normals.push([0.0, 0.0, 1.0]);
              }
              let normal = this.normal_map[[0.0, 0.0, 1.0]];
              this.faces.push([firstTop, neighborFirstTop, secondTop, normal]);
              this.faces.push([neighborFirstTop, firstTop, neighborSecondTop, normal]);
    
              normal = this.computeNormal([this.points[secondTop], this.points[neighborFirstTop]]);
              this.faces.push([secondTop, neighborFirstTop, secondTop + this.points.length, normal]);
              this.faces.push([secondTop + this.points.length, neighborFirstTop, neighborFirstTop + this.points.length, normal]);
    
              normal = this.computeNormal([this.points[neighborSecondTop], this.points[firstTop]]);
              this.faces.push([firstTop, firstTop + this.points.length, neighborSecondTop + this.points.length, normal]);
              this.faces.push([firstTop, neighborSecondTop + this.points.length, neighborSecondTop, normal]);
            }
          }
        }
      }
    
    NodeProcess() {
        for (let node of Object.values(this.nodes)) {
          if (node.neighbors.length == 1) {
            let first = node.points[0];
            let second = node.points[1];
            let normal = this.computeNormal([this.points[second], this.points[first]]);
            this.faces.push([first, second, first + this.points.length, normal]);
            this.faces.push([first + this.points.length, second, second + this.points.length, normal]);
          } else if (node.neighbors.length > 2) {
            let c = node.points[0];
            if (!([0.0, 0.0, 1.0] in this.normal_map)) {
              this.normal_map[[0.0, 0.0, 1.0]] = this.normals.length;
              this.normals.push([0.0, 0.0, 1.0]);
            }
            let normal = this.normal_map[[0.0, 0.0, 1.0]];
            for (let i = 1; i < node.points.length - 1; i++) {
                this.faces.push([c, node.points[i], node.points[i + 1], normal]);
            }
          }
        }
    }
    get() {
        let result = "";
        for (let p of this.points) {
            result += `v ${p[0]} ${p[1]} ${this.height}\n`;
          }
          for (let p of this.points) {
            result += `v ${p[0]} ${p[1]} 0\n`;
          }
          for (let vn of this.normals) {
            result += `vn ${vn[0]} ${vn[1]} ${vn[2]}\n`;
          }
          for (let f of this.faces) {
            result += `f ${f[0] + 1 + this.offset}//${f[3] + 1 + this.normal_offset} ${f[1] + 1 + this.offset}//${f[3] + 1 + this.normal_offset} ${f[2] + 1 + this.offset}//${f[3] + 1 + this.normal_offset}\n`;
          }
          return result;
    }
}

class Plain {
    constructor(data, base_index, offset, normal_offset) {
        this.height = data['height'];
        this.points = [];
        this.contour = [];
        this.normals = [];
        this.normal_map = {};
        this.point_index = {};
        this.segments = [];
        this.holes = [];
        this.faces = [];
        this.base_index = base_index;
        this.offset = offset;
        this.normal_offset = normal_offset;
        this.len_nodes = data['vertex'].length;
        this.loadData(data);
    }
    
    loadData(data) {
        let index = 0;
        let next = this.base_index;
        let hm = new Set();
        while (!hm.has(next)) {
            hm.add(next);
            const p = [data['vertex'][next - this.base_index]['x'], data['vertex'][next - this.base_index]['y']];
            this.point_index[p] = this.points.length;
            this.points.push(p);
            this.contour.push(new poly2tri.Point(p[0],p[1]));
            this.segments.push([index, (index + 1) % data['vertex'].length]);
            index += 1;
            if (!hm.has(data['vertex'][next - this.base_index]['neighbor'][0])) {
                next = data['vertex'][next - this.base_index]['neighbor'][0];
            } else {
                next = data['vertex'][next - this.base_index]['neighbor'][1];
            }
        }
        for (let hollow of data['hollow']) {
            let base = this.points.length;
            let hole = [];
            for (let i = 0; i < hollow.length; i++) {
                this.point_index[[hollow[i][0], hollow[i][1]]] = this.points.length;
                this.points.push([hollow[i][0], hollow[i][1]]);
                this.segments.push([base + i, base + (i + 1) % hollow.length]);
                hole.push(new poly2tri.Point(hollow[i][0], hollow[i][1]));
            }
            base += hollow.length;
            this.holes.push(hole);
        }
    }

    computeNormal(points) {
        const x1 = points[0][0];
        const y1 = points[0][1];
        const x2 = points[1][0];
        const y2 = points[1][1];
        const vector = [x2 - x1, y2 - y1];
        let normal = [-vector[1], vector[0]];
        const magnitude = Math.sqrt(normal[0] ** 2 + normal[1] ** 2);
        normal = [normal[0] / magnitude, normal[1] / magnitude];
        const key = [normal[0].toFixed(4), normal[1].toFixed(4), 0.0];
        if (!this.normal_map[key]) {
            this.normal_map[key] = this.normals.length;
            this.normals.push(key);
        }
        return this.normal_map[key];
    }

    getDirect() {
        let sum = 0;
        for (let i = 0; i < this.segments.length; i++) {
            let x1 = this.points[this.segments[i][0]][0];
            let y1 = this.points[this.segments[i][0]][1];
            let x2 = this.points[this.segments[i][1]][0];
            let y2 = this.points[this.segments[i][1]][1];
            sum += (x1 * y2 - x2 * y1);
        }
        return sum;
    }

    triangulateTop() {
        if (!([0.0, 0.0, 1.0] in this.normal_map)) {
            this.normal_map[[0.0, 0.0, 1.0]] = this.normals.length;
            this.normals.push([0.0, 0.0, 1.0]);
        }
        let normal = this.normal_map[[0.0, 0.0, 1.0]];
        let swctx = new poly2tri.SweepContext(this.contour);
        let triangles = swctx.addHoles(this.holes).triangulate().getTriangles();
        const self = this;
        triangles.forEach(function(t) {
            const f = [self.point_index[[t.getPoint(0).x, t.getPoint(0).y]], self.point_index[[t.getPoint(1).x, t.getPoint(1).y]], self.point_index[[t.getPoint(2).x, t.getPoint(2).y]], normal];
            self.faces.push(f);
        });
    }
        
    reconstruct() {
        this.triangulateTop();
        if (this.height !== 0) {
            let dir = this.getDirect();
            for (let i = 0; i < this.segments.length; i++) {
                let segment = this.segments[i];
                let normal;
                if (dir * this.height > 0) {
                    normal = this.computeNormal([this.points[segment[1]], this.points[segment[0]]]);
                } else {
                    normal = this.computeNormal([this.points[segment[0]], this.points[segment[1]]]);
                }
                if (dir > 0) {
                    this.faces.push([segment[0], segment[0] + this.points.length, segment[1], normal]);
                    this.faces.push([segment[0] + this.points.length, segment[1] + this.points.length, segment[1], normal]);
                } else {
                    this.faces.push([segment[0], segment[1], segment[0] + this.points.length, normal]);
                    this.faces.push([segment[0] + this.points.length, segment[1], segment[1] + this.points.length, normal]);
                }
            }
        }
        return [this.len_nodes, (this.height === 0) ? this.points.length : 2 * this.points.length, this.normals.length];
    }

    get() {
        let result = "";
        if (this.height != 0) {
            for (let p of this.points) {
                result += `v ${p[0].toFixed(4)} ${p[1].toFixed(4)} ${this.height.toFixed(4)}\n`
            }
        }
        for (let p of this.points) {
            result += (`v ${p[0].toFixed(4)} ${p[1].toFixed(4)} 0\n`);
        }
        for (let vn of this.normals) {
            result += `vn ${vn[0]} ${vn[1]} ${vn[2]}\n`;
        }
        for (let f of this.faces) {
            result += (
                `f ${
                f[0] + 1 + this.offset
                }//${f[3] + 1 + this.normal_offset} ${f[1] + 1 + this.offset}//${
                f[3] + 1 + this.normal_offset
                } ${f[2] + 1 + this.offset}//${f[3] + 1 + this.normal_offset}\n`
            );
        }
        return result;
    }
}

function ANSdecode(buffer, quant_bits = 12, renorm_bits = 16, window_bit = 16) {
    const alphabet = [];
    const freqs = [];
  
    const ft = new DataView(buffer.slice(0, 512));
    const codesReader = new DataView(buffer.slice(512));

    for (let i = 0; i < 256; i++) {
      const freq = ft.getUint16(i * 2, true);
      if (freq > 0) {
        alphabet.push(i);
        freqs.push(freq);
      }
    }

    const cdf = [];
    const mask = (1 << quant_bits) - 1;
    for (let i = 0; i < freqs.length + 1; i++) {
      cdf.push(freqs.slice(0, i).reduce((a, b) => a + b, 0));
    }

    let bufferLen = codesReader.byteLength;
    let outBuf = [];
    let code = codesReader.getUint32(0, true);
    bufferLen -= 4;

    while (code >= 1 << (2 * renorm_bits - window_bit)) {
        const s = code & mask;
        const index = cdf.findIndex((cdfi) => cdfi > s) - 1;
        code = freqs[index] * (code >>> quant_bits) + (code & mask) - cdf[index];
        outBuf.push(alphabet[index]);
        while (code < 1 << (2 * renorm_bits - window_bit) && bufferLen > 0) {
          code = (code << window_bit) + codesReader.getUint16(codesReader.byteLength - bufferLen, true);
          code = code>>>0;
          bufferLen -= 2;
        }
      }
    
      return new Uint8Array(outBuf.reverse());
}

function decode(binary_data) {
    let result = [];
    let decode_result = {'data': result, 'obj_data': []};
    let index = 0;
    let dv = new DataView(binary_data.buffer);
    let v_max = dv.getFloat32(index, false);
    index += 4;
    let v_min = dv.getFloat32(index, false);
    index += 4;
    while (index < binary_data.byteLength) {
      let header = dv.getUint16(index, false);
      index += 2;
      let data_type = header >> 14;
      let data_len = header & 0x3fff;
      if (data_type === 0) {
        let offset = dv.getUint16(index, false);
        index += 2;
        let height = dv.getFloat32(index, false);
        index += 4;
        let thickness = dv.getFloat32(index, false);
        index += 4;
        let vertex = [];
        for (let i = 0; i < data_len; i++) {
          let x = dv.getUint16(index, false);
          index += 2;
          let y = dv.getUint16(index, false);
          index += 2;
          let neighbor = [];
          while (true) {
            let n = dv.getUint16(index, false);
            index += 2;
            if (n === 0) {
              break;
            }
            neighbor.push(n + offset - 1);
          }
          vertex.push({'x': x, 'y': y, 'neighbor': neighbor});
        }
        buildEdge(vertex, offset);
        result.push({
          'type': 0,
          'height': height,
          'thickness': thickness,
          'vertex': vertex
        });
      } else if (data_type === 1) {
        let offset = dv.getUint16(index, false);
        index += 2;
        let height = dv.getFloat32(index, false);
        index += 4;
        let vertex = [];
        for (let i = 0; i < data_len; i++) {
          let x = dv.getUint16(index, false);
          index += 2;
          let y = dv.getUint16(index, false);
          index += 2;
          let neighbor = [];
          while (true) {
            let n = dv.getUint16(index, false);
            index += 2;
            if (n === 0) {
              break;
            }
            neighbor.push(n + offset - 1);
          }
          vertex.push({'x': x, 'y': y, 'neighbor': neighbor});
        }
        let hollow = [];
        while (true) {
          let hollow_len = dv.getUint8(index);
          index += 1;
          if (hollow_len === 0) {
            break;
          }
          let hollow_points = [];
          for (let j = 0; j < hollow_len; j++) {
            let x = dv.getUint16(index, false);
            index += 2;
            let y = dv.getUint16(index, false);
            index += 2;
            hollow_points.push([x, y]);
          }
          hollow.push(hollow_points);
        }
        buildEdge(vertex, offset);
        result.push({'type': data_type, 'height': height, 'vertex': vertex, 'hollow': hollow});
      } else if (data_type === 2) {
        result.push({'type': data_type, 'vertex_len': data_len});
      } else if (data_type === 3) {
        let vertex = [];
        for (let i = 0; i < data_len; i++) {
            let iindex = dv.getUint8(index);
            index += 1;
            let x = dv.getUint16(index, false);
            index += 2;
            let y = dv.getUint16(index, false);
            index += 2;
            let rotation = dv.getUint16(index, false);
            index += 2;
            let scale = dv.getFloat32(index, false);
            index += 4;
            vertex.push({'x': x, 'y': y, 'index': iindex, 'rotation': rotation, 'scale': scale});
        }
        let obj_len = dv.getUint8(index);
        index += 1;
        let obj_data = [];
            for (let i = 0; i < obj_len; i++) {
            let obj_offset = dv.getUint32(index, false);
            index += 4;
            let obj_length = dv.getUint32(index, false);
            index += 4;
            obj_data.push({'offset': obj_offset, 'length': obj_length});
        }
        result.push({'type': data_type, 'vertex': vertex, 'data': obj_data});
        while (index < dv.byteLength) {
            decode_result['obj_data'].push(dv.getUint8(index));
            index += 1;
        }
      }
    }
    dequantify(result, v_max, v_min);
    return decode_result;
}

function dequantify(data, v_max, v_min) {
    const factor = (v_max - v_min) / 0xffff;
    for (const d of data) {
      if (d.type !== 2) {
        for (const vertex of d.vertex) {
          vertex.x = vertex.x * factor + v_min;
          vertex.y = vertex.y * factor + v_min;
        }
      }
      if (d.type === 1) {
        for (const hollow of d.hollow) {
          for (const point of hollow) {
            point[0] = point[0] * factor + v_min;
            point[1] = point[1] * factor + v_min;
          }
        }
      }
      if (d.type === 3) {
        for (const vertex of d.vertex) {
          vertex.rotation = vertex.rotation * (360 / 0xffff);
        }
      }
    }
}
  
function buildEdge(vertexs, offset) {
    for (let i = 0; i < vertexs.length; i++) {
        for (const n of vertexs[i].neighbor) {
        if (n - offset < i) {
            vertexs[n - offset].neighbor.push(i + offset);
        }
        }
    }
}
  

async function reconstruction(filename, batch=true) {
    try {
        // const response = await fetch('decode.json');
        // const data = await response.json();
        const dracoLoader = new THREE.DRACOLoader();
        dracoLoader.setDecoderConfig({type:'js', decodeSpeed: 3, quantization: [ 11, 8, 8, 8, 8 ]});
        const response = await fetch(filename);
        const blob = await response.blob();
        const buffer = await new Promise(resolve => {
            const reader = new FileReader();
            reader.onload = () => resolve(reader.result);
            reader.readAsArrayBuffer(blob);
        });
        const data = decode(ANSdecode(buffer));
        let offset = 0;
        let base_index = 0;
        let normal_offset = 0;
        let structure = "";
        let objs = {obj: [], base_obj: []};
        let obj_index = [];
        if (batch) {
          for (let d of data['data']) {
              if (d['type'] === 0) {
                  const m = new Wall(d, base_index, offset, normal_offset);
                  const [len_nodes, len_points, len_normal] = m.reconstruct();
                  structure += m.get();
                  base_index += len_nodes;
                  offset += len_points;
                  normal_offset += len_normal;
              } else if (d['type'] === 1) {
                  const m = new Plain(d, base_index, offset, normal_offset);
                  const [len_nodes, len_points, len_normal] = m.reconstruct();
                  structure += m.get();
                  base_index += len_nodes;
                  offset += len_points;
                  normal_offset += len_normal;
              } else if (d['type'] === 2) {
                  if (d['vertex_len']) {
                      base_index += d['vertex_len'];
                  } else if (d['vertex']) {
                      base_index += d['vertex'].length;
                  }
              } else if (d['type'] === 3) {
                  objs.obj = d['vertex'];
                  obj_index = d['data'];
              }
          }          
        } else {
          for (let d of data['data']) {
              if (d['type'] === 0) {
                  const m = new Wall(d, base_index, 0, 0);
                  const [len_nodes, len_points, len_normal] = m.reconstruct();
                  structure += m.get();
                  structure += ";"
                  base_index += len_nodes;
                  offset += len_points;
                  normal_offset += len_normal;
              } else if (d['type'] === 1) {
                  const m = new Plain(d, base_index, 0, 0);
                  const [len_nodes, len_points, len_normal] = m.reconstruct();
                  structure += m.get();
                  structure += ";"
                  base_index += len_nodes;
                  offset += len_points;
                  normal_offset += len_normal;
              } else if (d['type'] === 2) {
                  if (d['vertex_len']) {
                      base_index += d['vertex_len'];
                  } else if (d['vertex']) {
                      base_index += d['vertex'].length;
                  }
              } else if (d['type'] === 3) {
                  objs.obj = d['vertex'];
                  obj_index = d['data'];
              }
          }
        }

        for (let i = 0; i < obj_index.length; i++) {
            const offset = obj_index[i].offset;
            const length = obj_index[i].length;
            const buffer = new Uint8Array(data['obj_data'].slice(offset, offset + length)).buffer;
            await dracoLoader.loadbuffer(buffer, function (geo) {
                objs.base_obj.push(geo);
            });
        }
        return {"structure": structure, "objs": objs};
    } catch (error) {
        console.log(error);
        return "";
    }
}
