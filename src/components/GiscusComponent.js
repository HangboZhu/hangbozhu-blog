import React from 'react';
import Giscus from "@giscus/react";
import { useColorMode } from '@docusaurus/theme-common';

export default function GiscusComponent() {
  const { colorMode } = useColorMode();
  const light_url = 'https://cdn.jsdelivr.net/gh/wjwei-handsome/CDN@0.3/comment.css';
  const dark_url = 'https://cdn.jsdelivr.net/gh/wjwei-handsome/CDN@0.3/comment-dark.css';
  const theme_url = colorMode == 'dark' ? dark_url : light_url;
  return (
    <div style={{marginTop:'30px'}}>
    <Giscus
      repo="HangboZhu/hangbozhu-blog"
      repoId="R_kgDOL9Vk4Q"
      category="General"
      categoryId="DIC_kwDOL9Vk4c4Cfdbp"
      mapping="pathname"
      term="Welcome to @giscus/react component!"
      strict="0"
      reactionsEnabled="1"
      emitMetadata="0"
      inputPosition="bottom"
      theme={theme_url} // 弄两套css，判断colormode.
      lang="en"
      loading="lazy"
      crossorigin="anonymous"
      async
    />
    </div>
  );}
