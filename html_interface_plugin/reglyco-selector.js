const tmpl = document.createElement('template');

tmpl.innerHTML = `
<style>
:host {
    width: 100%;
    height: 100%;
}
reglyco-select {
    overflow-y: auto;
    width: 100%;
    height: 30em;
    --max-select-display: 19;
    --palette-closer-style: var(--palette-closer-light);
}
reglyco-select option {
  display: none;
}

#residue_select {
    display: grid;
    grid-template-columns: repeat(10, 1fr);
    grid-template-rows: repeat(auto-fill, 1fr);
    max-height: 300px;
    width: 100%;
    overflow-x: hidden;
    overflow-y: auto;
}
#residue_select label {
    cursor: pointer;
    display: flex;
    flex-direction: row;
    font-size: 8pt;
    font-family: Helvetica, arial, sans-serif;
    border-radius: 5px;
    border: solid rgba(0,0,0,0) 1px;
}

#residue_select label:hover {
    border-color: #aaa;
}
#residue_select label[selected] {
    border-color: black;
}
#residue_select label input {
    position: absolute;
    top: -1000px;
    left: -1000px;
}
#residue_select label reglyco-sviewer {
    max-width: 50px;
    max-height: 50px;
    min-height: 50px;
    width: 100%;
    height: 100%;
    left: 0px;
}

section {
    display: block;
}
</style>

<section>
  <form id="residue_select">
  </form>
  <reglyco-select links id="sugar_select"></reglyco-select>
</section>
`;


const tmpl_input = document.createElement('template');

tmpl_input.innerHTML = `
  <label><reglyco-sviewer links renderer="svg"></reglyco-sviewer><span></span><input type="radio" name="residue"></input></label>
`;

const tmpl_sugar = document.createElement('template');

tmpl_sugar.innerHTML = `
  <option></option>
`;


let rewrite_seq_input = (inseq) => {
    let sequence = inseq;
    if (sequence.match(/5([AG])c/)) {
        sequence = sequence.replace(/5([AG])c/,'$1c');
        return rewrite_seq_input(sequence)
    }

    if (sequence.match(/[\)\]]([^[\)\]]+)([2346])S/)) {
        sequence = sequence.replace(/([\)\]])([^[\)\]]+)([2346])S/,'$1[HSO3(u?-$3)]$2');
        return rewrite_seq_input(sequence)
    }
    if (sequence.match(/([^[\)\]]+)([2346])S/)) {
        sequence = sequence.replace(/([^[\)\]]+)([2346])S/,'HSO3(u?-$2)$1');
        return rewrite_seq_input(sequence)
    }
    if (sequence.match(/([^[\)\]]+)([2346])Me/)) {
        sequence = sequence.replace(/([^[\)\]]+)([2346])Me/,'[Me(u?-$2)]$1');
        return rewrite_seq_input(sequence)
    }

    if (inseq !== sequence) {
        // console.log(inseq,sequence);
    }
    return sequence
};

let rewrite_seq_output = (outseq) => {
    if (! outseq) {
        return outseq;
    }
    let sequence = outseq;
    while (sequence.match(/Neu[AG]c/)) {
        sequence = sequence.replace(/Neu([AG])c/,'Neu5$1c')
    }
    while (sequence.match(/\[?HSO3\(u.-(\d)\)\]?([A-Za-z0-9]+)/)) {
        sequence = sequence.replace(/\[?HSO3\(u.-(\d)\)\]?([A-Za-z0-9]+)/,'$2$1S');
    }
    while (sequence.match(/\[?Me\(u.-(\d)\)\]?([A-Za-z0-9]+)/)) {
        sequence = sequence.replace(/\[?Me\(u.-(\d)\)\]?([A-Za-z0-9]+)/,'$2$1Me');
    }

    if (sequence.match(/[2346]S/)) {
        console.log(outseq,sequence);
    } else {
        if (outseq.indexOf('HSO3') >= 0) {
            console.log(outseq,sequence);
        }
    }

    if (outseq !== sequence) {
        // console.log(outseq,sequence);
    }
    return sequence
};


class MySugarClass extends window.SViewer.IupacSugar {
    set sequence(inputseq) {
        super.sequence = rewrite_seq_input(inputseq);
    }
    get sequence() {
        return rewrite_seq_output(super.sequence);
    }
}

class MySelect extends window.SugarSelect.default {
    connectedCallback() {
        this.SugarClass = MySugarClass;
        super.connectedCallback();
    }
}

class MyViewer extends window.SViewer.SViewerLite {
    connectedCallback() {
        this.SugarClass = MySugarClass;
        super.connectedCallback();
    }
}


let values = {};

function wire_events() {
  this.shadowRoot.querySelector('#residue_select').addEventListener('change', (ev) => {
    if ( ! ev.target || ev.target.getAttribute('name') !== 'residue' ) {
        return;
    }
    this.updateSugars();
    this.shadowRoot.querySelector('reglyco-select').refreshOptions();
    let identifier = this.shadowRoot.querySelector('#residue_select').residue.value;
    let [residue,chain,aa] = identifier.split('_');
    if (values[identifier]) {
        this.shadowRoot.querySelector('reglyco-select').value = values[identifier];
    } else {
        this.shadowRoot.querySelector('reglyco-select').reset();
    }
    if (this.shadowRoot.querySelector('[selected]')) {
        this.shadowRoot.querySelector('[selected]').removeAttribute('selected');
    }
    ev.target.parentNode.setAttribute('selected','');
  });
  this.shadowRoot.querySelector('reglyco-select').addEventListener('change', async () => {
    let identifier = this.shadowRoot.querySelector('#residue_select').residue.value;
    let [residue,chain,aa] = identifier.split('_');
    let seq = this.shadowRoot.querySelector('reglyco-select').value
    values[identifier] = seq;
    let label = this.shadowRoot.querySelector(`#residue_select input[value='${identifier}']`).parentNode;
    label.querySelector('reglyco-sviewer').SugarClass = MySugarClass;
    label.querySelector('reglyco-sviewer').sequence = seq;
    label.querySelector('reglyco-sviewer').scaleToFit();
  });
}

class ReGlycoSelector extends HTMLElement {

  #config = {}

  async connectedCallback() {
    let shadowRoot = this.attachShadow({mode: 'open'});
    shadowRoot.appendChild(tmpl.content.cloneNode(true));
    wire_events.call(this);
  }

  set configuration(config) {
    this.#config = config;
    this.reconfigure();
    let temp_selector = document.createElement('reglyco-select');
    for (let aa of ['SER','THR','ASN','TRP']) {
        let sugars = this.#config.configurations[aa].map( glycan => glycan.iupac );
        MySelect.processSequences(sugars,MySugarClass);
    }
  }

  get value() {
    return values;
  }

  updateSugars() {
    let [residue,chain,aa] = this.shadowRoot.querySelector('#residue_select').residue.value.split('_');
    let sugars = this.#config.configurations[aa].map( glycan => glycan.iupac );
    this.shadowRoot.querySelector('#sugar_select').innerHTML = '';
    for (let sugar of sugars) {
      let kid = tmpl_sugar.content.cloneNode(true);
      kid.firstElementChild.innerText = sugar;
      this.shadowRoot.querySelector('#sugar_select').appendChild(kid);
    }
  }

  glytoucanMap() {
    return Object.fromEntries(Object.entries(this.#config.configurations).map( ([aa,glycans]) => glycans ).flat().map( glycan => [glycan.iupac, glycan.glytoucan] ));
  }

  reconfigure() {
    this.shadowRoot.querySelector('#residue_select').innerHTML = '';
    for (let info of this.#config.glycosylation.available) {
      let residue = info.residueID;
      let chain = info.residueChain;
      let aa = info.residueName;
      let kid = tmpl_input.content.cloneNode(true);
      kid.querySelector('input').setAttribute('value',`${residue}_${chain}_${aa}`);
      kid.querySelector('input').setAttribute('residue',`${residue}`);
      kid.querySelector('input').setAttribute('chain',`${chain}`);
      kid.querySelector('input').setAttribute('aa',`${aa}`);
      kid.querySelector('input').innerText = `${residue}_${chain}`;
      kid.querySelector('label span').appendChild(document.createTextNode(`${aa} ${residue}`));
      kid.querySelector('label').style.gridRow = ( Math.floor(parseInt(residue) / 10) ) + 1;
      kid.querySelector('label').style.gridColumn = ( parseInt(residue) % 10 ) + 1;
      this.shadowRoot.querySelector('#residue_select').appendChild(kid);
    }
  }
}

customElements.define('reglyco-sviewer',MyViewer);
customElements.define('reglyco-select',MySelect);
customElements.define('reglyco-selector',ReGlycoSelector);

