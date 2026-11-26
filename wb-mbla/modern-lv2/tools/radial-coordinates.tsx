// Convert a hollow-fibre wall distance to the complementary radius and fixed CP units.
export function RadialCoordinates({ y, ratio, exact }: { y: number; ratio?: number; exact?: boolean }) {
  const normalized = (value: number) => exact ? String(value) : Number(value.toPrecision(6)).toString();
  return <>
    <span data-y-over-r>{`y/R = ${normalized(y)}`}</span> /{' '}
    <span data-r-over-r>r/R = {Number((1 - y).toPrecision(8))}</span>
    {ratio !== undefined && <> / <span data-y-over-cp>y/δ<sub>CP</sub> = {Number((y / ratio).toPrecision(8))}</span></>}
  </>;
}
