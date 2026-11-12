/* Minimal SVG recording context for the spatial renderer. No external dependencies. */
(function(root) {
  'use strict';
  const escape = value => String(value).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&apos;'}[c]));
  const num = value => { if (!Number.isFinite(value)) throw new Error('Invalid vector coordinate'); return String(+value.toFixed(5)); };
  class SVGContext {
    constructor(minFont=11) {
      this.items=[];this.defs=[];this.stack=[];this.path='';this.clips=[];this.minFont=minFont;
      this.fillStyle='#000';this.strokeStyle='#000';this.lineWidth=1;this.globalAlpha=1;
      this.font='12px "Times New Roman"';this.textAlign='left';this.lineCap='butt';this.dash=[];
    }
    setTransform() {} // Draw commands use logical coordinates, not device pixels.
    clearRect() {}
    save() { this.stack.push(Object.fromEntries(['fillStyle','strokeStyle','lineWidth','globalAlpha','font','textAlign','lineCap','dash','clips'].map(k=>[k,Array.isArray(this[k])?[...this[k]]:this[k]]))); }
    restore() { Object.assign(this,this.stack.pop()); }
    beginPath() { this.path=''; }
    moveTo(x,y) { this.path+=`M${num(x)} ${num(y)} `; }
    lineTo(x,y) { this.path+=`L${num(x)} ${num(y)} `; }
    closePath() { this.path+='Z '; }
    rect(x,y,w,h) { this.moveTo(x,y);this.lineTo(x+w,y);this.lineTo(x+w,y+h);this.lineTo(x,y+h);this.closePath(); }
    arc(x,y,r,start,end) {
      if (Math.abs(end-start-Math.PI*2)>1e-6) throw new Error('Only full-circle spatial markers are supported');
      this.path+=`M${num(x+r)} ${num(y)} A${num(r)} ${num(r)} 0 1 0 ${num(x-r)} ${num(y)} A${num(r)} ${num(r)} 0 1 0 ${num(x+r)} ${num(y)} Z `;
    }
    clip() { const id='clip'+this.defs.length;this.defs.push(`<clipPath id="${id}"><path d="${this.path}"/></clipPath>`);this.clips.push(id); }
    add(item) { for(const id of this.clips)item=`<g clip-path="url(#${id})">${item}</g>`;this.items.push(item); }
    paint(stroke) {
      this.add(`<path d="${this.path}" fill="${stroke?'none':escape(this.fillStyle)}" opacity="${num(this.globalAlpha)}"`+
        (stroke?` stroke="${escape(this.strokeStyle)}" stroke-width="${num(this.lineWidth)}" stroke-linecap="${escape(this.lineCap)}" stroke-dasharray="${this.dash.map(num).join(' ')}"`:'')+'/>');
    }
    fill() { this.paint(false); }
    stroke() { this.paint(true); }
    fillRect(x,y,w,h) { this.add(`<rect x="${num(x)}" y="${num(y)}" width="${num(w)}" height="${num(h)}" fill="${escape(this.fillStyle)}" opacity="${num(this.globalAlpha)}"/>`); }
    strokeRect(x,y,w,h) { const saved=this.path;this.beginPath();this.rect(x,y,w,h);this.stroke();this.path=saved; }
    setLineDash(dash) { this.dash=[...dash]; }
    fillText(text,x,y) {
      const size=Math.max(this.minFont,parseFloat(this.font)||12);
      this.add(`<text x="${num(x)}" y="${num(y)}" font-family="Times New Roman, Times, serif" font-size="${num(size)}" fill="${escape(this.fillStyle)}" text-anchor="${this.textAlign==='center'?'middle':this.textAlign==='right'?'end':'start'}">${escape(text)}</text>`);
    }
    drawImage(image,x,y,w,h) {
      let uri=image.src;
      if(!uri.startsWith('data:image/')){
        const copy=document.createElement('canvas');copy.width=image.naturalWidth;copy.height=image.naturalHeight;
        copy.getContext('2d').drawImage(image,0,0);uri=copy.toDataURL('image/png');
      }
      this.add(`<image x="${num(x)}" y="${num(y)}" width="${num(w)}" height="${num(h)}" preserveAspectRatio="none" opacity="${num(this.globalAlpha)}" href="${uri}"/>`);
    }
    markup() { return `<defs>${this.defs.join('')}</defs>`+this.items.join(''); }
  }
  root.LUCENT_VECTOR={SVGContext,escape,num};
  if(typeof module==='object'&&module.exports)module.exports=root.LUCENT_VECTOR;
})(typeof window!=='undefined'?window:globalThis);
