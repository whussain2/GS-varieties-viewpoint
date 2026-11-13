import sys, asyncio
from playwright.async_api import async_playwright
src, out = sys.argv[1], sys.argv[2]   # html path, output stem
async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch()
        pg = await b.new_page(viewport={"width": 680, "height": 1000}, device_scale_factor=4)
        await pg.goto("file://" + src); await pg.wait_for_timeout(600)
        h = await pg.evaluate("document.querySelector('.page').getBoundingClientRect().height")
        await pg.set_viewport_size({"width": 680, "height": int(h) + 2})
        await pg.screenshot(path=out + ".png", full_page=True)
        await pg.pdf(path=out + ".pdf", width="180mm", height=f"{(h+2)*25.4/96:.1f}mm", print_background=True)
        print("height px", h)
        await b.close()
asyncio.run(main())
