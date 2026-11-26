// Adapt the current profile selection to the research site's own component.
import { createRoot } from 'react-dom/client';
import { FiltrationSchematic } from './schematic';
import type { RadialScale } from './coordinates';

declare global {
  interface Window {
    wbMountSiteSelector: (
      element: HTMLElement,
      onSelect: (z: number, wallDistance: number) => void,
    ) => (z: number, wallDistance: number, radialScale: RadialScale) => void;
  }
}

window.wbMountSiteSelector = (element, onSelect) => {
  const root = createRoot(element);
  return (z, wallDistance, radialScale) => {
    root.render(<FiltrationSchematic geometry="HF" position={z}
      wallDistance={wallDistance} radialScale={radialScale} onSelect={onSelect} />);
  };
};
